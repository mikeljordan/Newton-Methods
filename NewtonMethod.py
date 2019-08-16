#Newton's Method for solving non-linear system
from numpy import *
import pylab as pl
def f(x):   #define functions f1=f[0], f2=f[1]
	f=zeros((len(x)),dtype=float64)
	f[0]=log(x[0]**2+x[1]**2)-sin(x[0]*x[1])-log(2)-log(pi)
	f[1]=exp(x[0]-x[1]) + cos(x[0]*x[1])
	return f  #return values of f at x
def jacobian1(f,x): #compute Jacobian
	jac1=array([[((2*x[0])/(x[0]**2+x[1]**2))-x[1]*cos(x[1]*x[0]),
	((2*x[1])/(x[0]**2+x[1]**2))-x[0]*cos(x[1]*x[0])],
	[exp(x[0]-x[1])-x[1]*sin(x[0]*x[1]),-exp(x[0]-x[1])-x[0]*sin(x[0]*x[1])]])
	f0=f(x)
	return jac1,f0  #return the value of jacobian and value of f at x
en = [] # list of error (en)
def newtonRaphson3(f,x,tol=1.0e-6):
	print "Approximated solutions"	
	for i in range(1,30):
		jac1,f0 = jacobian1(f,x)
#solve the non-linear system 'jac1(roots)=-f0"		
		roots = linalg.solve(jac1,-f0) 
		x = x + roots				
		error = linalg.norm(roots,inf) #compute the error using inf.norm
		en.append(error)
		print x #list the approximated solutions
		if (error) < tol:return x,i
if __name__ == "__main__":
	en1 = [] # list of error (en+1)
	x = array([2., 2.]) #initial guess
	roots,i = newtonRaphson3(f,x)
	print "Roots=", roots
	print "Iterations=",i
	for i in range(len(en)-1):
		en1.append(en[i+1])
	en.remove(en[-1])
	en = pl.log(en) #compute the log of errors in the list en
	en1 = pl.log(en1) #compute the log of errors in the list en+1
	print "Error list(en):\n",en 
	print "Error list(en+1):\n",en1
	pl.plot(en,en1)
	slope, intercept = polyfit(en,en1, 1)
	print "Speed of convergence(slope of Graph)=",slope 
	print "Intercept=",intercept
	pl.title("Graph of log(en)=plog(en1)+c")
	pl.show()
