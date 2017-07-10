
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np

def Symbol_Maker(X) :

	"""
	
	Generates an array of different discrete colours for plotting. 
	
	(Perhaps need more of a continuous generalization) 

	"Extended Colmaker does both colours and symbols??"
	
	Parameters: 
	
	X: array
	
	
	Returns:
	
	Col_Ret: array of floats
	
	
	
	"""

	Symbols = [ "o" , "v" , "s", "*" , "h" , "+" , "<"]
	Cols = [ 'r' , 'b' , 'g' , 'm' , 'c' , 'k' , 'y' ]
	#Cols = [ 'r' ] 
	Col_Ret = [ ]
	
	"""
	Temprorarily remove colours and symbols:
	for i in range( len(X) ) :

		
		K = (i+1) % len(Symbols)
		for j in range( 7 ) :
		    J = (j+1) % 7
		    Col_Ret.append( Cols[ 0 ] + Symbols[K] ) 
	
	"""	 
	
	for i in range( len(X) ) :
		K = (i+1) % len(Symbols)
		Col_Ret.append( Cols[K] + Symbols[K] ) 
		
	return Col_Ret
	
def Colour_Maker(X) :

	"""
	
	Generates an array of different discrete colours for plotting. 
	
	(Perhaps need more of a continuous generalization) 

	"Extended Colmaker does both colours and symbols??"
	
	Parameters: 
	
	X: array
	
	
	Returns:
	
	Col_Ret: array of floats
	
	
	
	"""

	Symbols = [ "o" , "v" , "s", "*" , "h" , "+" , "<"]
	Cols = [ 'r' , 'b' , 'g' , 'm' , 'c' , 'k' , 'y' ]
	#Cols = [ 'r' ] 
	Col_Ret = [ ]
	
	"""
	Temprorarily remove colours and symbols:
	for i in range( len(X) ) :

		
		K = (i+1) % len(Symbols)
		for j in range( 7 ) :
		    J = (j+1) % 7
		    Col_Ret.append( Cols[ 0 ] + Symbols[K] ) 
	
	"""	 
	
	for i in range( len(X) ) :
		K = (i+1) % len(Cols)
		Col_Ret.append( Cols[K]  ) 
		
	return Col_Ret

def Histogram(Sample, filename , bin_num , x_axis , y_axis = 'Frequency' , Show_Plot = False , Colour = 'r') :

		#Set up the axis text sizes:
		
		matplotlib.rc('xtick', labelsize=20) 
		matplotlib.rc('ytick', labelsize=20) 
		
		
		#Set up the fond size:
		
		font = {
        'weight' : 'bold',
        'size'   : 18}
		
		#This code makes sure the axis labels stay within the plot:
		from matplotlib import rcParams
		rcParams.update({'figure.autolayout': True})
		plt.rc('font', **font)
		plt.rc('axes', labelsize = 28 )
		#Make the figure: 
		plt.figure(1)
		#plt.tight_layout()
		plt.clf()
		fig1 = open(filename + '.png' , 'w')
		n, bins, patches = plt.hist(Sample, bin_num,  normed=1,alpha=0.75, histtype='stepfilled')
		plt.setp(patches, 'facecolor', Colour, 'alpha', 0.75)
		Half_Bin_width = (bins[2] - bins[1])/2
		bins = [x+Half_Bin_width for x in bins]
		plt.hold(True)
		plt.plot(bins[0:(len(bins)-1)],n,'k')
		#plt.title('Algebraic Connectivity Distribution for Conditioned Networks')
		plt.xlabel( x_axis )
		plt.ylabel(y_axis) 
			
		plt.savefig(fig1 , format ='png' )
		fig1.close()
		if Show_Plot == True :
			plt.show()



def Compare_Histogram(Sample, Sample_2 , filename , bin_num , x_axis , y_axis = 'Frequency' , Show_Plot = False , Colour = 'r') :

		#Set up the axis text sizes:
		matplotlib.rc('xtick', labelsize=17) 
		matplotlib.rc('ytick', labelsize=17) 
		plt.rc('axes', labelsize = 20 )
		
		
		n, bins, patches = plt.hist(Sample, bin_num,  normed=1,alpha=0.75, histtype='stepfilled')
		plt.setp(patches, 'facecolor', Colour, 'alpha', 0.75)
		Half_Bin_width = (bins[2] - bins[1])/2
		bins = [x+Half_Bin_width for x in bins]
		
		n2, bins2, patches2 = plt.hist(Sample_2, bin_num,  normed=1,alpha=0.75, histtype='stepfilled')
		plt.setp(patches2, 'facecolor', Colour, 'alpha', 0.75)
		Half_Bin_width2 = (bins2[2] - bins2[1])/2
		bins2 = [x+Half_Bin_width2 for x in bins2]
		
		
		
		
		
		#plt.hold(True)
		plt.figure(1)
		plt.clf()
		fig1 = open(filename + '.png' , 'w')
		plt.plot(bins[0:(len(bins)-1)],n,'k--' , linewidth=4)
		
		
		plt.hold(True)
		
		plt.plot(bins2[0:(len(bins2)-1)],n2,'g' , linewidth=4)
		
		#plt.title('Algebraic Connectivity Distribution for Conditioned Networks')
		plt.xlabel( x_axis )
		plt.ylabel(y_axis) 
			
		plt.savefig(fig1 , format ='png' )
		fig1.close()
		if Show_Plot == True :
			plt.show()

def Plot_Multiple( X_Data , Y_Data , Y_Errors = 0 , With_Error_Bars = True ) :

	return 0


def Error_Plot(x_variables,y_variables, y_errs , filename , colours ,x_lab , y_lab, Legend_Text , LogX = False , LogY = False, Legend_Box_Below = False , Show_Plot = False , Axis = None ) :
	#Take in x and y variables as [ ] 's containing several different data sets to be plotted
	
	
	#Set up the axis text sizes:
		
	matplotlib.rc('xtick', labelsize=20) 
	matplotlib.rc('ytick', labelsize=20) 


	#Set up the fond size:
	"""
	font = {
    'weight' : 'bold',
    'size'   : 18}
	"""
	#This code makes sure the axis labels stay within the plot:
	from matplotlib import rcParams
	rcParams.update({'figure.autolayout': True})
	#plt.rc('font', **font)
	plt.rc('axes', labelsize = 28 )
	
	
	
	
	
	
	#Number of X and Y variables 
	X_plots = len(x_variables)
	#Open file
	plt.figure(1)
	plt.clf()
	fig2 = open(filename + '.png' , 'w' )
	
	Figure_List = [ ]
	#Plotting loop
	for i in range(X_plots) : 
		plt.plot( x_variables[i] , y_variables[i] , colours[i] )
		Current_Handle = plt.errorbar( x_variables[i] , y_variables[i] , y_errs[i] , fmt = colours[i] + 'o' )
		Figure_List.append( Current_Handle ) 
		
	plt.xlabel( x_lab )
	plt.ylabel( y_lab ) 
	if LogX == True :
		plt.xscale('log')
		
	if LogY == True : 
		plt.yscale('log') 
	
	print("Making Plot")
	plt.tight_layout()
	#Save Figure
	
	if Axis != None : 
		plt.axis( Axis ) 
	
	
	if Legend_Box_Below == True :
		LEG = plt.legend(Figure_List, Legend_Text , loc='upper center', bbox_to_anchor=(0.5, -0.1),
		  fancybox=True, shadow=True , ncol=3)
		plt.savefig(fig2 , format = 'png' , bbox_extra_artists=(LEG,), bbox_inches='tight')
	else :
		plt.legend(Figure_List, Legend_Text , loc = 1 )
		#Save Figure
		plt.savefig(fig2 , format = 'png' )

	if Show_Plot == True :
		plt.show()

	fig2.close()	
	

def Normal_Plot(x_variables,y_variables , filename , colours ,x_lab , y_lab, Legend_Text , LogX = False , LogY = False , Legend_Box_Below = False ,Show_Plot = False ) :
	#Take in x and y variables as [ ] 's containing several different data sets to be plotted
	
	
	#Set up the axis text sizes:
		
	matplotlib.rc('xtick', labelsize=20) 
	matplotlib.rc('ytick', labelsize=20) 
	
	
	#Set up the fond size:
	
	font = {
    'weight' : 'bold',
    'size'   : 18}
	
	#This code makes sure the axis labels stay within the plot:
	from matplotlib import rcParams
	rcParams.update({'figure.autolayout': True})
	plt.rc('font', **font)
	plt.rc('axes', labelsize = 28 )
	
	
	#Number of X and Y variables 
	X_plots = len(x_variables)
	
	#Open file
	plt.figure(1)
	plt.clf()
	fig2 = open(filename + '.png' , 'w' )
	
	Figure_List = [ ]
	
	#Plotting loop
	for i in range(X_plots) : 
		plt.plot( x_variables[i] , y_variables[i] , colours[i] )
		Current_Handle, = plt.plot( x_variables[i] , y_variables[i]  , colours[i] + 'o' )
		Figure_List.append( Current_Handle ) 
	
		
	
	if LogX == True :
		plt.xscale('log')
		
	if LogY == True : 
		plt.yscale('log') 
		
	plt.xlabel( x_lab )
	plt.ylabel( y_lab ) 
	
	if Legend_Box_Below == True :
		LEG = plt.legend(Figure_List, Legend_Text , loc='upper center', bbox_to_anchor=(0.5, -0.1),
		  fancybox=True, shadow=True , ncol=3)
		plt.savefig(fig2 , format = 'png' , bbox_extra_artists=(LEG,), bbox_inches='tight')
	else :
		plt.legend(Figure_List, Legend_Text , loc = 1 )
		#Save Figure
		plt.savefig(fig2 , format = 'png' )
		
	if Show_Plot == True :
			plt.show()
		
	fig2.close()
	
def Scatter_Plot(x_variables,y_variables , filename , colours ,x_lab , y_lab, Legend_Text , LogX = False , LogY = False , Show_Plot = False , Legend_Box_Below = False) :
	#Take in x and y variables as [ ] 's containing several different data sets to be plotted
	
	
	
	#Number of X and Y variables 
	X_plots = len(x_variables)
	
	#Open file
	plt.figure(1)
	plt.clf()
	fig2 = open(filename + '.png' , 'w' )
	
	Figure_List = [ ]
	
	#Plotting loop
	for i in range(X_plots) : 
		Current_Handle, = plt.plot( x_variables[i] , y_variables[i]  , colours[i] + 'o' )
		Figure_List.append( Current_Handle ) 
	
	plt.legend(Figure_List, Legend_Text , loc = 1 )	
		
	plt.xlabel( x_lab )
	plt.ylabel( y_lab ) 
	if LogX == True :
		plt.xscale('log')
	if LogY == True : 
		plt.yscale('log') 
	
	
	
	if Legend_Box_Below == True :
		LEG = plt.legend(Figure_List, Legend_Text , loc='upper center', bbox_to_anchor=(0.5, -0.1),
		  fancybox=True, shadow=True , ncol=3)
		plt.savefig(fig2 , format = 'png' , bbox_extra_artists=(LEG,), bbox_inches='tight')
	else :
		plt.legend(Figure_List, Legend_Text , loc = 1 )
		#Save Figure
		plt.savefig(fig2 , format = 'png' )
	fig2.close()
	
	if Show_Plot == True :
		plt.show()

def Plot_Heat_Map(X_arr,Y_arr, M , filename , x_lab , y_lab , ColBar_Label = ' '  )  : 

	#Set up the text of the axis:

	X, Y = np.meshgrid(X_arr, Y_arr)

	plt.figure(3)
	
	"""
	#plt.figure(num=None, figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k' ) 
	
	fig = open( filename + '.png'  , 'w' ) 
	plt.clf()
	#fig2, ax = plt.subplots()
	
	plt.pcolor(X,Y,np.swapaxes(M,0,1))
	#plt.contour(X,Y,np.swapaxes(M,0,1), cmap=plt.cm.bone)
	ColBar = plt.colorbar()
	#plt.clim(0,0.5)
	#Set the scale of the colour bar so that it is the same for both periodic and solid boundary condition cases.
	ColBar.set_label(ColBar_Label)
	ax = plt.gca()
	ax.set_xticks(np.arange(M.shape[1]) + 0.5, minor=False)
	ax.set_yticks(np.arange(M.shape[0]) + 0.5, minor=False)
	
	#ax.set_xticks(np.arange(M.shape[1]) , minor=False)
	#ax.set_yticks(np.arange(M.shape[0]) , minor=False)
	"""
	
	fig, ax = plt.subplots()
	#fig.subplots_adjust(bottom=0.25,left=0.25) # make room for labels
	plt.rc('axes', labelsize = 28 )
	heatmap = ax.pcolor(X,Y,np.swapaxes(M,0,1))
	cbar = plt.colorbar(heatmap)

	# Set ticks in center of cells
	#Adjusted from: https://stackoverflow.com/questions/32980633/adding-text-ticklabels-to-pcolor-heatmap 
	ax.set_xticks(np.arange(M.shape[0]) + 0.5, minor=False)
	ax.set_yticks(np.arange(M.shape[1]) + 0.5, minor=False)

	# Rotate the xlabels. Set both x and y labels to headers[1:]
	ax.set_xticklabels(Y_arr[1:])
	ax.set_yticklabels(X_arr[1:])

	
	plt.xlabel(x_lab)
	plt.ylabel(y_lab)
	plt.savefig( filename, format = 'png' ) 
	#fig.close()
	
	
def Plot_Properties_ColourBar( Property_1 , Property_2 , Colour_Property, filename , xlabel , ylabel , col_label) : 
	plt.figure()
	plt.clf()
	fig2 = open(filename + '.png' , 'w' )
	#X_Array , Y_Array, Z_Array, Legend_Text = self.Sort_Via_Property( Property_1 , Property_2 , Colour_Property ) 
	X_Array = Property_1
	Y_Array = Property_2
	Z_Array = Colour_Property
	#Col = Col_Maker(Z_Array) 
	cm = plt.cm.get_cmap('RdYlBu')
	#print( X_Array )
	#print( len(X_Array ) ) 
	sc = plt.scatter( X_Array, Y_Array, c = Z_Array , cmap = cm ) 
	colbar = plt.colorbar(sc)
	colbar.set_label( col_label ) 
	plt.xlabel(xlabel)
	plt.ylabel( ylabel ) 
	plt.savefig(fig2 , format = 'png' )
	#plt.show()
	fig2.close()

"""

Function : dobuleplot - plot two things in the same box??



"""




#Do we need a sperate plotting library for network plotting functions (and for example stuff that relies on plotly??? ) 
