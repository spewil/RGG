
import numpy as np
import itertools
#Kappa_Array = [ 20.0 ]
Kappa_Array = np.arange(25.0 , 50.0 , 5.0 ) 
N_array = [ 1000 ]
d_array = np.arange(2,11,1)
#d_array = [ 5 ] 
BC_Array = [ 'S' ]


Samples_Per_Submission = 100
Submissions = 1


Que_Type = 'standard' 

#Open the master submit file and go to the correct file:
Master_Submit_File = open( "Master_Submit_RGG.sh" , 'a' ) 
Master_Submit_File.write('cd /scratchcomp03/MGarrod_4\n')

for i,j,k,l in itertools.product(Kappa_Array, N_array,d_array,BC_Array):
	#This is very bad code. We keep swapping the order of N and Kappa. Tragic. 
	filename = "Script_Submit_Kap_" + str(i) + "_N_" + str(j) + "_d_" + str(k) + "_BC_" + str(l)
	Script = open( filename + '.sh' , 'a' ) 
	Script.write( "#!/bin/bash\n#PBS -N RGG_" + str(j) + "_" + str(i) + "_" + str(k) + "_" + str(l) + "\n#PBS -q standard\ncd /scratchcomp03/MGarrod_4/\n" + 
	"chmod u+x /scratchcomp03/MGarrod_4/RGG_Sample_Mu2_Only.py\npython /scratchcomp03/MGarrod_4/RGG_Sample_Mu2_Only.py " + str(j) + " " + str(i) + " " +str(k) + " " + str(l) + " " + str(Samples_Per_Submission) ) 
	Script.close()
	
	
	#Before submitting scripts we must write a command to move into the working directory:
	#(NB: could generalize this code to a function which specifies the working directory)
	
	for p in range( Submissions ) : 
		Master_Submit_File.write( "qsub " + filename + '.sh\n' ) 

#Close files:		
Master_Submit_File.close()	

#Now create master submission script to submit everything
#Remember we would like to be able to do things in parallel!

