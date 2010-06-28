import sys
for i in sys.argv:
	print i
infile = sys.argv[1]
outfile = sys.argv[2]
try:
	f = open(infile,"r")
except IOError:                     
	print "The file does not exist, exiting gracefully"
xRef=0
yRef=0
actualX=0
actualRow=[]
positionMatrix=[]
counter =0
actualSegment=""
level=False
print "this is the file"+str(infile)
for line in f:
	if line.startswith("tDescription"):
		description=line.strip().split(" ")
		size=description[5]
		size=size.split("x")
		width=size[0]
		height=size[1]
		row=["size",str(width),str(height)]
		positionMatrix.append(row)
	if line.startswith("[Level0]"):
		level=True
	if line.startswith("[") and level==True:
		row=[]
		splittedLine1=line.strip().split("[")
		splittedLine2=splittedLine1[1].strip().split("]")
		name=splittedLine2[0]
	if line.startswith("x") and level==True:
		splittedLine1=line.strip().split("=")
		xValue=splittedLine1[1]
	if line.startswith("y") and level==True:
		splittedLine1=line.strip().split("=")
		yValue=splittedLine1[1]
		row=[name,str(xValue),str(yValue)]
		positionMatrix.append(row)
	
f.close()
outfile = open(outfile,"w")

for row in positionMatrix:
	outfile.write("\t".join(row)+"\n")
outfile.close()
