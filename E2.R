install.packages("seqinr")
library(seqinr)
help(seqinr)
aedes = read.fasta(file="aedesaegypti.fasta")
aedes = aedes[[1]]  
n = length(aedes)

################### Question 1 ################################

################### a)
m = 200
k = n%/%m
gcc = numeric(k)
for(i in 1:k){
  a = (i-1)*m+1
  b = a+m
  gcc[i] = GC(aedes[a:b])
}
ts.plot(gcc)

m = 1000
k = n%/%m
gcc2 = numeric(k)
for(i in 1:k){
  a = (i-1)*m+1
  b = a+m
  gcc2[i] = GC(aedes[a:b])
}
ts.plot(gcc2)

################### b)
pos_max_gcc = which.max(gcc)
max_gcc = gcc[pos_max_gcc]
print(pos_max_gcc)
print(max_gcc)
ts.plot(gcc[70611:72611])

pos_max_gcc2 = which.max(gcc2)
max_gcc2 = gcc2[pos_max_gcc2]
print(pos_max_gcc2)
print(max_gcc2)
ts.plot(gcc[259435:261435])

################### Question 2 ################################

################### Estimate the probabilities
prob = table(aedes)/length(aedes)
prob
# count(aedes,1,freq=T)

################### Transition probability matrix
a1=count(aedes[1:77706756],2)
a2=count(aedes[77706757:155413512],2)
a3=count(aedes[155413513:233120268],2)
a4=count(aedes[233120269:310827022],2)
a = a1+a2+a3+a4
c = matrix(a, 4, 4, byrow=TRUE, dimnames = list(c("a", "c", "g", "t"),c("a", "c", "g", "t")))
c[,]/(c[,1]+c[,2]+c[,3]+c[,4])

################### Multinomial model 
n=length(aedes)
par=3
c1=count(aedes[1:103609007],1)
c2=count(aedes[103609008:207218015],1)
c3=count(aedes[207218016:310827022],1)
c = c1+c2+c3
p1=count(aedes[1:103609007],1,freq=T)
p2=count(aedes[103609008:207218015],1,freq=T)
p3=count(aedes[207218016:310827022],1,freq=T)
p = p1+p2+p3
BIC=-2*sum(c*log(p)) + par*log(n)
BIC # 162.642.198 Multinomial is better!!!

################### Markov chain model 
n=length(aedes)-1
par=12
a1=count(aedes[1:77706756],2)
a2=count(aedes[77706757:155413512],2)
a3=count(aedes[155413513:233120268],2)
a4=count(aedes[233120269:310827022],2)
a = a1+a2+a3+a4
c = matrix(a, 4, 4, byrow=TRUE, dimnames = list(c("A","C","G","T"),c("A","C","G","T")))
p=c[,]/(c[,1]+c[,2]+c[,3]+c[,4])
BIC=-2*sum(c*log(p)) + par*log(n)
BIC #838.882.695

################### Question 3 ################################

# Is there any relationship between the presence of "acc" and the GC?
n=length(aedes)
m = 200
k = n%/%m
gcc = numeric(k)
acc = numeric(k)
ex3 = numeric(k)
for(i in 1:k){
  a = (i-1)*m+1
  b = a+m-1
  gcc[i] = GC(aedes[a:b])
  acc[i] = count(aedes[a:b],3)["acc"]
  if (acc[i]!=0){
    ex3[i]=1
  }
  else{
    ex3[i]=0
  }
}
ts.plot(gcc)
ts.plot(ex3)

# What is the probability of "acc" for a chunk with a GC content of 0.51? 
# 1st chunk
y1 = ex3[1:388534]
x1 = gcc[1:388534]
fit = glm(y1 ~ x1, family=binomial)
summary(fit)
e1= coef (summary(fit))[, "Estimate"]
s1= coef (summary(fit))[, "Std. Error"]

# 2nd chunk
y1 = ex3[388535:777068]
x1 = gcc[388535:777068]
fit = glm(y1 ~ x1, family=binomial)
e2= coef (summary(fit))[, "Estimate"]
s2= coef (summary(fit))[, "Std. Error"]

# 3rd chunk
y1 = ex3[777069:1165602]
x1 = gcc[777069:1165602]
fit = glm(y1 ~ x1, family=binomial)
e3= coef (summary(fit))[, "Estimate"]
s3= coef (summary(fit))[, "Std. Error"]

# 4th chunk
y1 = ex3[1165603:1554135]
x1 = gcc[1165603:1554135]
fit = glm(y1 ~ x1, family=binomial)
e4= coef (summary(fit))[, "Estimate"]
s4= coef (summary(fit))[, "Std. Error"]

sw=1/(s1^2) + 1/(s2^2) + 1/(s3^2) + 1/(s4^2)
est1=(e1/(s1^2) + e2/(s2^2) + e3/(s3^2) + e4/(s4^2))/sw
se1=1/(sw^.5)
est1
se1

Y1 = est1[1] + 0.51*est1[2]
Y1
P1 = plogis(Y1)
P1

# Plot the estimated probability of "acc" against the GC content.
Y11 = numeric(100)
P11 = numeric(100)
gcc_plot = seq(0, 0.99, by=0.01)
for (i in (1:100)) {
  Y11[i]=est1[1] + gcc_plot[i]*est1[2]
  P11[i]=plogis(Y11[i])
}
Y11
P11
plot(gcc_plot,P11)

################### Question 4 ################################

# Is there any relationship between the mean of the counts of "aaa" and the GC content? 
n=length(aedes)
m = 200
k = n%/%m
gcc = numeric(k)
aaa = numeric(k)
for(i in 1:k){
  a = (i-1)*m+1
  b = a+m-1
  gcc[i] = GC(aedes[a:b])
  aaa[i] = count(aedes[a:b],3)["aaa"]
}
ts.plot(gcc)
ts.plot(aaa)

# What is the predicted mean of the counts of "aaa" for a chunk with a GC content of 0.4 ?
# 1st chunk
y2 = aaa[1:388534]
x2 = gcc[1:388534]
fit2 = lm(y2 ~ x2)
summary(fit2)
e21= coef (summary(fit2))[, "Estimate"]
s21= coef (summary(fit2))[, "Std. Error"]

# 2nd chunk
y2 = aaa[388535:777068]
x2 = gcc[388535:777068]
fit2 = lm(y2 ~ x2)
e22= coef (summary(fit2))[, "Estimate"]
s22= coef (summary(fit2))[, "Std. Error"]

# 3rd chunk
y2 = aaa[777069:1165602]
x2 = gcc[777069:1165602]
fit2 = lm(y2 ~ x2)
e23= coef (summary(fit2))[, "Estimate"]
s23= coef (summary(fit2))[, "Std. Error"]

# 4th chunk
y2 = aaa[1165603:1554135]
x2 = gcc[1165603:1554135]
fit2 = lm(y2 ~ x2)
e24= coef (summary(fit2))[, "Estimate"]
s24= coef (summary(fit2))[, "Std. Error"]

sw2=1/(s21^2) + 1/(s22^2) + 1/(s23^2) + 1/(s24^2)
est2=(e21/(s21^2) + e22/(s22^2) + e23/(s23^2) + e24/(s24^2))/sw2
se2=1/(sw2^.5)
est2
se2

Y2 = est2[1] + 0.4*est2[2]
Y2

# Plot the estimated mean of the counts of "aaa" against the GC content.
Y22 = numeric(100)
gcc_plot = seq(0, 0.99, by=0.01)
for (i in (1:100)) {
  Y22[i]=est2[1] + gcc_plot[i]*est2[2]
}
Y22
plot(gcc_plot,Y22)

