GAUSS NEWTON SOLVER
-------------------

THERE ARE 4 PROGRAMS CONTAINED IN THIS ZIP FILE FOR PART 1

gn_exp.java
gn_log.java
gn_qua.java
gn_rat.java

Each of these programs have the same way of executing them.

These are Command Line programs so you will first have to open up Terminal or Command Prompt if you're using Windows. Once you have this opened up, using the 'cd' statement, move to the current directory of where you extracted the zip file. For example:

	cd home/user/Desktop/Calc3Project/GaussNewton

		or if using Windows:

	cd user\Desktop\Calc3Project\GaussNewton

Now that you are in this current directory type:

	java -cp .:Jama-1.0.3.jar gn_exp

		or if you are using Windows type:

	java -cp .;Jama-1.0.3.jar gn_exp

This should run the program.
Now you can call the same line but change the ending to gn_qua or gn_rat, depending on which program
you want to run.


1)
	When you run one of the programs you will be prompted with a question like this:

		This is the Gauss-Newton solver for a Quadratic function.
		Please type in the name of the text file containing a list of points. (ex. data.txt)

	This is where you will need to write in the path of the text file you are trying it to read. For
	example:

		/home/user/Desktop/data.txt

		or if you are using Windows:

		C:\joe\Desktop\data.txt



2)
	After this you will be prompted with:

		What is your initial guess for parameter a?

		What is your initial guess for parameter b?

		What is your initial guess for parameter c?

	Here is where you can put in your initial guesses each time being prompted



3)
	Proceeding this you will be prompted with:

		How many iterations do you want to run the Gauss-Newton algorithm?

	and you can type in the number of iterations you want.



4)
	The final thing the program will prompt you for is:

		Which QR Factorization method would you like to run?
		Householder or Givens. Type "H" or "G"

	Here is where you can type either 'H' or 'G', and it is not case-sensitive. It will either choose
	to calculate QR by Householder Reflections or Givens rotations.

	It will also give you the Q and R matricies for each iteration.


	If you choose 'G' for instance the ouput you will recieve looks something like this:

		Using the Givens method, this converges to:
		0.16
	 	-2.0
	 	0.85

----------------------------------------------------------------------------------------------------------



	



	










