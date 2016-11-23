import os, re, sys, yaml
from subprocess import call
from optparse import OptionParser

################################################################################################################################


################################################################################################################################
# Parameters
################################################################################################################################

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)


parser.add_option("-i","--inputfile",
                  dest = "inputfile",
                  default = "",
                  metavar = "path",
                  help = "Input file [required]")

parser.add_option("-o","--outputfile",
                  dest = "outputfile",
                  default = "",
                  metavar = "path",
                  help = "Name of the output file (by default it is the same as the input) [required]")

parser.add_option("-t","--template",
                  dest = "template",
                  default = "",
                  metavar = "path",
                  help = "Pandoc template file to use (if not specified use pandoc default) [required]")

parser.add_option("-T","--title",
                  dest = "title",
                  default = "Untitled tutorial",
                  metavar = "string",
                  help = "Title of tutorial in quotes (if not already in yaml header) [required]")

parser.add_option("-S","--subtitle",
                  dest = "subtitle",
                  default = "",
                  metavar = "string",
                  help = "Subtitle of tutorial in quotes (if not already in yaml header) [required]")

parser.add_option("-r","--removefiles",
				  action = "store_false",
				  default = True,
                  dest = "remove",              
                  metavar = "boolean",
                  help = "Remove temporary files (default=%default) [required]")

parser.add_option("-L","--latex",
				  action = "store_true",
				  default = False,
                  dest = "latex",              
                  metavar = "boolean",
                  help = "Compile output file with pdflatex (default=%default) [required]")


(options,args) = parser.parse_args()

inputfile  = os.path.abspath(options.inputfile)
outputfile = os.path.abspath(options.outputfile) if options.outputfile != "" else inputfile[:inputfile.rfind('.')]+".tex"
template   = os.path.abspath(options.template)

title      = options.title
subtitle   = options.subtitle

remove     = options.remove
latex      = options.latex
bibtex     = False


################################################################################################################################


def getPandocCall(inputfile, outputfile, template=""):
	call = ["pandoc"]

	# Output file
	call.append("-o")
	call.append(outputfile)

	# Standalone file
	call.append("-s")
	if (template != ""):
		call.append("--template")
		call.append(template)

	# Do not parse html blocks
	call.append("-R")

	# Use listings package for code blocks
	call.append("--listings")

	# Input file
	call.append(inputfile)

	return(call)
#


def getYamlHeader(text):

	match = re.match("(\A---.*?)---",text, re.DOTALL)

	header = yaml.load(match.groups(0)[0])
	if (title != "Untitled tutorial"):
		header["title"] = title

	if (subtitle != ""):
		header["subtitle"] = subtitle

	if ("author" in header.keys()):
		authors = header["author"].split(",")
		if (len(authors) > 1):
			header["author"] = ", ".join(authors[:-1]) + " and " + authors[-1]
	else:
		header["author"] = "Anonymous"

	return((header,text[match.end():]))
#


def parseLiquid(text, header=None):

	start = 0
	while (True):

		# Match the next liquid tag
		# {% (tag) (content) %}
		# Can be multiline
		match = re.search(r'{%\s(.*?)\s(.*?)%}', text, re.DOTALL)

		# No more matches
		if (match == None):
			break

		(tag,content) = match.groups()

		# Process tags
		if (tag == "eq"):
			replacement = "\\begin{equation}\n\t%s\n\end{equation}" % content.strip()
		elif (tag == "eqinline"):
			replacement = "`$ %s $`" % content.strip()
			#replacement = ""
		elif (tag == "cite"):
			replacement = "\citep{%s}" % content[:content.find("--")].strip()
		elif (tag == "bibliography"):			
			parts = content.split()
			for i in range(0,len(parts)):
				if (parts[i] == "--file"):
					bibfile = parts[i+1]

			if (bibfile.find('/') > 0):
				bibfile = bibfile[bibfile.rfind('/')+1:]

			if (header != None):
				header["bibtex"] = bibfile
				bibtex = True

			replacement = ""
		else:
			sys.stdout.write("WARNING Unsupported tag: %s\n" % (text[match.start():match.end()]))
			sys.stdout.write("skipping...\n")
			replacement = ""

		text = text[:match.start()] + replacement + text[match.end():]

	return(text)
#


def parseFigures(text):

	start = 0
	while (True):

		# Match the next html figure tag
		# Can be multiline
		match = re.search(r'<figure>.*?<a id="(.*?)".*?<img(.*?)>.*?<figcaption>(.*?)</figcaption>.*?</figure>', text, re.DOTALL)

		# No more matches
		if (match == None):
			break

		(label, figure, caption) = match.groups()

		# Process \includegraphics
		scale = "[width=\\textwidth]"
		for part in figure.strip().split():
			(tag,content) = part.split("=")

			if (tag == "src"):
				figfile = content.replace('"',"")

			if (tag == "style"):
				if (content[1:6].lower() == "width" and content[-3:-1] == '%;'):
					mult  = int(content[content.find(':')+1:content.find('%')])/100
					scale = "[width=%.6f\\textwidth]" % mult
				else:
					sys.stdout.write("WARNING Unsupported image style specification %s in '%.20s...'\n" % (part, caption))
					sys.stdout.write("skipping...\n")
		
		# Process \caption
		# Should actually run the caption through pandoc to parse the caption body text
		caption = caption.replace("%","\%")
		if (re.match(r'[f|F]igure\s\d+:',caption)):
			caption = caption[caption.find(':')+1:].strip()

		replacement  = "\\begin{figure}\n\t\centering\n"
		replacement += "\t\includegraphics%s{%s}\n" % (scale, figfile)
		replacement += "\t\caption{%s}\n" % caption.strip()
		replacement += "\t\label{%s}\n" % label.strip()
		replacement += "\end{figure}\n"

		text = text[:match.start()] + replacement + text[match.end():]

	return(text)
#


def parseFigureRefs(text):

	start = 0
	while (True):

		# Match the next figure reference
		# [Figure (number)](#(label))
		# Can be multiline
		match = re.search(r'\[Figure\s\d+\]\(#(.+?)\)', text, re.DOTALL)

		# No more matches
		if (match == None):
			break

		label = match.groups(0)[0]

		replacement = "Figure \\ref{%s}" % label.strip()
		text = text[:match.start()] + replacement + text[match.end():]

	return(text)
#


def removeMdLines(text):

	start = 0
	while (True):

		# Match the next horizontal line
		# Three or more dashes (---)
		# Can be multiline
		match = re.search(r'---*', text, re.DOTALL)

		# No more matches
		if (match == None):
			break

		text = text[:match.start()] + "\n\clearpage\n" + text[match.end():]

	return(text)
#


def formatInlineMath(text):

	start = 0
	while (True):

		# Match the next inline equation
		# \lstinline!$(math)$!	
		# Can be multiline
		match = re.search(r'\\lstinline\!\$(.+?)\$\!', text, re.DOTALL)

		# No more matches
		if (match == None):
			break

		label = match.groups(0)[0]

		replacement = "$ %s $" % label.replace("\\\\","\\").strip()
		text = text[:match.start()] + replacement + text[match.end():]

	return(text)
#


def formatSuperscript(text):

	super_regex = [r'(\\\^\{\}\((.+?)\))[^\)]', r'(\\\^\{\}(\w+?))\W']

	for regex in super_regex: 

		start = 0
		while (True):

			# Match the next text superscript (no parentheses)
			# \^{}(superscript)
			# Can be multiline
			#match = re.search(r'\\\^\{\}(\w+?)\W', text, re.DOTALL)
			#match = re.search(r'\\\^\{\}[\((.+?)\)[^\)]', text, re.DOTALL)
			match = re.search(regex, text, re.DOTALL)

			# No more matches
			if (match == None):
				break

			print(match.groups(1)[1])

			replacement = "$^{%s}$" % match.groups(1)[1]
			text = text[:match.start()] + replacement + text[match.end():]

	return(text)
#


def removeRefSection(text):

	start = 0
	while (True):

		# Remove the heading relevant references at the end
		# \section{Relevant References}\label{relevant-references}
		# Can be multiline
		match = re.search(r'\\section\{Relevant References\}\\label{relevant-references}', text, re.DOTALL & re.IGNORECASE)

		# No more matches
		if (match == None):
			break

		text = text[:match.start()] + text[match.end():]

	return(text)
#




################################################################################################################################

##################
# Pre-processing #
##################

text = open(inputfile,'r').read()

# Read header
(header, text) = getYamlHeader(text)

# Parse liquid tags
text = parseLiquid(text, header)

# Parse figures in html tags
text = parseFigures(text)

# Parse figure references
text = parseFigureRefs(text)

# Remove horizontal lines
text = removeMdLines(text)

# Add header to text and save temporary file
outfile = open(inputfile[:inputfile.rfind('.')]+"-temp"+inputfile[inputfile.rfind('.'):],'w', encoding='utf-8')
outfile.write("---\n")
#for key in header:
#	outfile.write("%s : %s\n" % (key, header[key]))
yaml.dump(header, outfile, default_flow_style=False, encoding='utf-8')
outfile.write("---\n")
outfile.write(text)
outfile.close()



#####################
# Pandoc conversion #
#####################
call(getPandocCall(inputfile[:inputfile.rfind('.')]+"-temp"+inputfile[inputfile.rfind('.'):], 
				   inputfile[:inputfile.rfind('.')]+"-temp.tex", template=template))




###################
# Post-processing #
###################

text = open(inputfile[:inputfile.rfind('.')]+"-temp.tex",'r').read()

# Framed boxes (quotes)
text = text.replace("\\begin{quote}","\\begin{framed}").replace("\end{quote}","\end{framed}")

# Format inline math (not properly converted by pandoc)
text = formatInlineMath(text)

# Format text superscripts (does not work perfectly)
text = formatSuperscript(text)

# Remove relevant references heading
text = removeRefSection(text) 

# Remove unknown latex tags that pandoc inserts
# remove \toprule, \bottomrule, \tightlist
text = text.replace("\\toprule","").replace("\\bottomrule","").replace("\\tightlist","")


outfile = open(outputfile,'w',encoding='utf-8')
outfile.write(text)
outfile.close()


# Remove temporary files




#################
# Compile latex #
#################

if (latex):
	call(['pdflatex',outputfile])
	if (bibtex):
		call(['bibtex',outputfile[:outputfile.rfind('.')]])		
		call(['pdflatex',outputfile])
	call(['pdflatex',outputfile])




