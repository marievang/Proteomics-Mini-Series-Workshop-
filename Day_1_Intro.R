#Proteomics 

#Day 1 - Introducing to R for Biostatistics

#R is extensible, easy to use, good for graphics, separates easy input data from operations

# R supports numeric variables as you encountered in algebra
a = 2
b = 5
c = a+b

# Enter the variable by itself to learn its value
c

# Not all variables are numeric, though.
# Boolean values are useful for making decisions

a = TRUE    #often abbreviated as a=T
b = FALSE

# logical AND: 
a&&b
# logical OR: 
a||b

#Numeric comparisons can produce Boolean results

5>=7
4==4      #We use "=" as the assignment operator.  "==" is a comparison!
4!=3


##################################################################################
#Data Types
#Vectors let us do many operations at once
#The c() function creates a vector

a = c(2,8,43,3,8)
b = c(21,7,9,34,4)
sums = a+b
sums                #adds together corresponding elements when we add vectors


#When vectors are different sizes, strange things can happen: is is fine as long as the longer object
#length is  not a multiple of shorter length
a = c(1,2,3,4,5,6)
b = c(0,1)
alternating = a+b
alternating


#We can pull out individual elements
sums[2] #takes the second element of sums vector.

sums[c(2,4)]

#Some operators and functions help us make vectors
sums[2:4] #all numbers from a to b (a:b)

evens = 2*(1:20)
evens

series = 14:20
series

repeats = rep(NA,10) #rep is repeating a specific value x times. 
repeats

randoms = runif(100) #uniform distribution, gives random samples of this distr.
randoms
hist(randoms)
#Many functions exist to help us look at vectors:

mean(randoms)      #arithmetic mean 
summary(randoms)   #more detailed description
length(randoms)    

plot(sort(randoms))

#Factors are very useful with categorical data, 

fruit = c("apple", "apple", "pear", "banana", "banana")
fruitf = factor(fruit)           

fruit
fruitf
#Without duplicating, print the distinct values ONLY 
levels(fruitf)



sex = c("m", "f", "m", "m", "f", "m")
sexf = factor(sex)
# By default, a factor displays as a barplot
plot(sexf)

#Matrices are two-dimensional tables of values (they must be the same type)
mat = matrix(c(1,2,3,4,5,6),byrow=T,nrow=2)
mat

#We can even name the rows and columns!
dimnames(mat) = list(c("row1","row2"), c("col1","col2","col3"))


#transpose
tmat = t(mat)

#append
row3 = c(7,8,9)
biggermat = rbind(mat,row3) #concat two things together 

#Distributions are available for us to generate random numbers
uniform = runif(100)
normal = rnorm(100)
poisson = rpois(100, lambda=5)

#We can envision these distributions through histograms
hist(uniform)
hist(normal)
hist(poisson)

#Many functions return lists as output: t.test example

a = rnorm(10, mean=3)
b = rnorm(10, mean=5)
c = t.test(a,b)

#What data type is 'c'?
typeof(c)

#What is the structure of 'c'?
str(c)

#The structure report tells us what we can extract from the t.test result: $ 
c$p.value

#Data frames are one of the most common types for people reading spreadsheets:

label = c(rep("case",5),rep("control",5))
reading = c(22, 25, 27, 28, 24, 12, 17, 13, 16, 17)
df = data.frame(label,reading)

df
### Functions
# If you are performing the same complex operation at many places within your code,
# you could instead define a function to carry it out.

# This function computes the distance between two points in Euclidean geometry
# The Pythagorean Theorem says that for a right triangle, A^2 + B^2 = C^2
# xdiff^2 + ydiff^2 = distance^2
Pythagoras = function(x1,y1, x2=0,y2=0) {
  xdiff = abs(x1-x2);
  ydiff = abs(y1-y2);
  distance = sqrt(xdiff^2 + ydiff^2);
  return(distance);
}

# The same function that finds the distance from the origin to a single point...
Pythagoras(3,4)

# can find the distance between two points
Pythagoras(1,6,4,10)


### Examine Transcript Data in Excel
# Drag df_transcripts.csv to an Excel window
# Walk through the column headers
# Show the CSV in a text editor
# Sort the Excel table by geneNames and marvel at Date interpretation for Septin-6.

### Read df_transcripts.csv into R
# Data published as https://doi.org/10.1038/nprot.2016.095
Transcripts <- read.table(file="df_transcripts.csv",sep=",",header=T)

# check if a file is in the working directory
#CHECK WHAT IS A P.VALUE AND WHAT IS A Q.VALUE

list.files()

# Functions for a helpful snapshot of data frames
summary(Transcripts)
head(Transcripts)
tail(Transcripts)

# Learn the size of a data frame
dim(Transcripts)

# Pull out a particular row of a data frame
# This line pulls out just row 13, grabbing all columns
Transcripts[13,]
# This line pulls out just rows 13-16, grabbing all columns
Transcripts[13:16,]

# Grab out rows matching a criterion (out of a big data object)
subset(Transcripts,geneIDs=="MSTRG.11")

# Pull out a particular column of a data frame
summary(Transcripts[18])
# Rather than referencing columns by number, we can use their names
summary(Transcripts$pval)

# Pull out a particular cell of a data frame
Transcripts[13,18]

# Pull out rows with p-value below 0.05
LowP=(subset(Transcripts,pval<0.05))
dim(LowP)

# Write these rows with low p-values to a new file
write.table(LowP,file="LowPValues.csv",sep=",",row.names=F)

# Add a new column, the log base 2 of fold change
log2fc = log2(Transcripts$fc)
ExpandedTranscripts = cbind(Transcripts,log2fc)
summary(ExpandedTranscripts)


### Examine data for a particular transcript yielding a low p-value
Target = Transcripts[191,]
Females = c(Target$FPKM.ERR188234,Target$FPKM.ERR188245,Target$FPKM.ERR188273,Target$FPKM.ERR188337,Target$FPKM.ERR188428,Target$FPKM.ERR204916)
Males =   c(Target$FPKM.ERR188044,Target$FPKM.ERR188104,Target$FPKM.ERR188257,Target$FPKM.ERR188383,Target$FPKM.ERR188401,Target$FPKM.ERR188454)

# Let's produce our own p-value, using the less-sophisticated T-test.
t.test(Males,Females)

# Plot the male and female FPKM values for this gene
#ylim: limits on the y axis
plot(Males,pch="m", ylim=c(0,110),xlab="person",ylab="FPKM") #pch: plotting character

points(Females,pch="f",col="red")

#comment on the data: 
#the highest value of the males is still lower than the lower value of the females.













