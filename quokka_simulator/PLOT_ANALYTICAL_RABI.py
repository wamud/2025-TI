import time
from datetime import datetime

print(datetime.now().strftime("%d/%m/%y %H:%M:%S\n"))

start = time.time()
debut = time.process_time()

import numpy as np
import matplotlib.pyplot as plt
import qutip
import sys
import math
from colorsys import hls_to_rgb
import re
import os
from collections import defaultdict

def giveMoreExactValues(θ,φ):
    # We will be subbing values of a=cos(θ/2) and b=e^(iφ)*sin(θ/2) (from the initial state of n rotated qubits (a|0> + b|1>)^n)
    # into the analytical equations for aL and bL in filea and fileb. Numpy introduces a lot of numerical errors though.
    # For example, sin(2*np.pi / 2) ≠ 0. Hence need to do a bit of tinkering here to get more exact values for a and b:

    #Here are some exact values for theta:

    if θ == 2*np.pi:
        cosθon2 = -1
        sinθon2 = 0
    elif θ == 0:
        cosθon2 = 1
        sinθon2 = 0
    elif θ == np.pi:
        cosθon2 = 0
        sinθon2 = 1
    elif θ == 0.5*np.pi:
        c = 0.7071067811865475244008443621048490392848359376884740365883398689953662392310535194251937671638207863675069231154561485124624180279253686063220607485499679157066113329637527963778999752505763910302857350547799858029851372672984310073642587093204445993047761646152421543571607254198813018139976257039948436266982731659044148203103076291761975273728751438799808649177876101687659285056771873017042494235801934499853495024075152720138951582271239115342464684593107902892315557983343565065078092844936186176442546324306247
        cosθon2 = c
        sinθon2 = c

    elif θ == 1.5*np.pi:
        c = 0.7071067811865475244008443621048490392848359376884740365883398689953662392310535194251937671638207863675069231154561485124624180279253686063220607485499679157066113329637527963778999752505763910302857350547799858029851372672984310073642587093204445993047761646152421543571607254198813018139976257039948436266982731659044148203103076291761975273728751438799808649177876101687659285056771873017042494235801934499853495024075152720138951582271239115342464684593107902892315557983343565065078092844936186176442546324306247
        cosθon2 = -c
        sinθon2 = c
    elif θ == 7*np.pi/8:
        cosθon2 = 0.195090322016128267848284868477022240927691617751954807754502089494763318785924580225325309234090381730992070105536611758965799834847607388826591434427552071609461968800368541406340990778672482876294780099080346776443730183027931917537989861258456667025731263924584964200153894836250465831928554723818708921835149190078236881898892695024212958043909589362065670501747286702576028923828164299618278957604139950704205164167518430815244864526180355930338185004302566840934726674186221019671544334876057352355431553003905672183128794774965050246082812166207693042719720630578951551354505204668896423033364176769180659222493609557976665298459281048029011718619960106013784519204975529213787027047295427053782651333205341776163755727803937852799490078182305998749998123551148530377469065069087472195521301492311040384816327009258028479523127416199097987040356088474065078658472684157740141479879686895104927241974061190750523979422127541540070611059889808850622298279215439993364661609020753712509076232693050412202553869705638482349490921098878100420463065242030222614186177355923069981280454077789198630442683704990349229824272896922786287636923014220438629817737159966172181890624476212971380311083163992137048984244149199680501065590706089650301025997669894486736316401134553832321619531553898439231121033760719027110288028415003241952499714718645314557134942046452587061814765272473422405635910416803702484472183715808316289682381693829662686515400835600835365617661549724387060029928841114340344022796107441002226162026092036830173987335795922302976517449581494430634823080901206820243337760648212969352538651475158035337737996038495474981126395491908856507480369138632771138786692570977604865083142256320937252527415311903841115560921824099383940638189283613814622259866012357660362355396501793938200661773535555606957212526218889038703563288291969736222836476342329568468413560763268239117013641737048566353584664646199191096703975860421735208054704771634922619973141304685234240305526439161200319363003855592731868734708239120719145253284710249936795160236938574132644137
        sinθon2 = 0.9807852804032304491261822361342390369739337308933360950029160885453065135496050639150649858533007632598948662798775784681310960848381701091485451909052981223580423918286860736338652741318972946739839332937486597435047390244869403252433116449612877676624274785105265658772764361560079035227698048820627724374873999548813877728693895193453570149596600565114596216808811598974974376637649736484339907745256748778850919666529636079121537087166579311692869625392730877273169182359525555967002080191631891777385320737152586760573313838758177972229360952418770644878194072961014142561469497983941173015281678513506977368248549117568681406296661672515780967848295595313183626978459266364655599277726571034191824647267066122286157397777810017310280906084874980384570416361485703534799986923778386115367396266286124759179569750440786872680071119173791762535549115348387964301228952946521613775493268121499714518794089237664232621492243933388271076478675294076659927687655752891775417643816463078769886597115145615435177901462367925444974705108848116707152307340310238739390557442102857543313360251970133052691748216064136906864520735657831681544433096797253795653256770007057090704374733196928541263641224763371566271421656871895916062786011162846545380849890444525135370632508538800300349197628448240221271159390231407146032768378316691612795668522124590960977442823530365569835296708891875755722950757618832430641346381978272284525288876796275570164418754763898634527764552359856187958505988157473827614802463618978944582588582967301208385707093222539875594256501722380264559567050690811024888510523327496337032954544999402033977754414176710618627463341200553732106701470713982314668003845606313544261067100621668803959665436814740182290546307513811481832355693826638447027578456422067372646679502054749857993518366982988908222022340405147341586081888969125095509043051366223547239750445994161571362716146060070274394680981444195465772489644992433114044920571858085365066818293679195239296860280529833397422606174648244104569986251289742797717764317056651198938748605307985488441762799620594753356231459683196574146707206989252094300785783538087750023979931429925821784749437212662315641801211254007627816205274788983641801842624165977271723767966983158688194256913775358497297023928310943965634980238603800048138728826592896178075639054029477302830755012110775065234933497986822343358627964500010176899458595377662685005375022370914533248705202435

    elif θ == np.pi/8:
        cosθon2 = 0.980785280403230449126182236134239036973933730893336095002916088545306513549605063915064985853300763259894866279877578468131096084838170109148545190905298122358042391828686073633865274131897294673983933293748659743504739024486940325243311644961287767662427478510526565877276436156007903522769804882062772437487399954881387772869389519345357014959660056511459621680881159897497437663764973648433990774525674877885091966652963607912153708716657931169286962539273087727316918235952555596700208019163189177738532073715258676057331383875817797222936095241877064487819407296101414256146949798394117301528167851350697736824854911756868140629666167251578096784829559531318362697845926636465559927772657103419182464726706612228615739777781001731028090608487498038457041636148570353479998692377838611536739626628612475917956975044078687268007111917379176253554911534838796430122895294652161377549326812149971451879408923766423262149224393338827107647867529407665992768765575289177541764381646307876988659711514561543517790146236792544497470510884811670715230734031023873939055744210285754331336025197013305269174821606413690686452
        sinθon2 = 0.195090322016128267848284868477022240927691617751954807754502089494763318785924580225325309234090381730992070105536611758965799834847607388826591434427552071609461968800368541406340990778672482876294780099080346776443730183027931917537989861258456667025731263924584964200153894836250465831928554723818708921835149190078236881898892695024212958043909589362065670501747286702576028923828164299618278957604139950704205164167518430815244864526180355930338185004302566840934726674186221019671544334876057352355431553003905672183128794774965050246082812166207693042719720630578951551354505204668896423033364176769180659222493609557976665298459281048029011718619960106013784519204975529213787027047295427053782651333205341776163755727803937852799490078182305998749998123551148530377469065069087472195521301492311040384816327009258028479523127416199097987040356088474065078658472684157740141479879686895104927241974061190750523979422127541540070611059889808850622298279215439993364661609020753712509076232693050412202553869705638482349490921098878100420463065242030222614186177355923069981280454077789198630442683704990349229824272896922786287636923014220438629817737159966172181890624476212971380311083163992137048984244149199680501065590706089650301025997669894486736316401134553832321619531553898439231121033760719027110288028415003241952499714718645314557134942046452587061814765272473422405635910416803702484472183715808316289682381693829662686515400835600835365617661549724387060029928841114340344022796107441002226162026092036830173987335795922302976517449581494430634823080901206820243337760648212969352538651475158
    else:
        cosθon2 = np.cos(θ/2)
        sinθon2 = np.sin(θ/2)
    
    #And some exact values for phi:


    if φ == np.pi:
        expterm = -1
    elif φ == 0:
        expterm = 1
    elif φ == 2*np.pi:
        expterm = 1
    elif φ == np.pi/2:
        expterm = 1j
    elif φ == 3*np.pi/2:
        expterm = -1j
    else:
        expterm = np.exp(1j*φ)

    a = cosθon2
    b = expterm * sinθon2

    return a,b

def calcLogicalCoords(αL,βL):
    # This function converts an αL and βL from αL|0〉+ βL|1〉to find 
    # θL and φL with the state in the form cos(θL/2)|0_L〉+ e^(iφL)sin(θL/2)|1_L〉
    # Then finds the spherical coords x,y,z to plot on qutip bloch sphere.
    if αL == βL == 0:
        print("Unreachable trajectory")
        return 0,0,0,0

    # Now convert this aL and bL to logical theta, logical phi and x,y,z:
    
    aLr = αL.real
    aLi = αL.imag

    r_alpha = math.sqrt((aLr**2) + (aLi**2))
    if r_alpha >  1: # Occasionally comes out as 1.0000000000000002 so round to 1
        r_alpha = 1

    theta = 2 * np.arccos(r_alpha) #THES ES LOGICAL THETA!!

    r_beta = math.sqrt((βL.real**2) + (βL.imag**2))
    theta2 = 2*np.arcsin(r_beta) # should be the same as theta, but just in case

    # Added checks just in case
    # print("thetaL calculated from alphaL ",theta/np.pi,"π")
    # print("thetaL calculated from betaL ",theta2/np.pi,"π")
    
    if theta > 0 and theta2 < 0:
        print("thetaL's have different signs")
    if theta < 0 and theta2 > 0:
        print("thetaL's have different signs")
    if abs(theta) - abs(theta2) > 1e-6:
        print("θL calculated from αL ≠ θL calculated from βL. \n θLα = ",theta/np.pi,"π, θLβ = ",theta2/np.pi,"π")


    #Manually set some values to avoid numerical error.
    #If logical theta is 0, for example, then the logical coordinate is just (x,y,z)=(0,0,1).

    if theta == 0 or theta == 2*np.pi:
        x = 0
        y = 0
        z = 1
        phi = 0

    elif theta == np.pi:
        x = 0
        y = 0
        z = -1
        phi = 0
    
    else:

        #Now to figure out logical phi, via each θ from aL and bL.
        # First convert the imaginary numbers aL and bL to their polar coords
        # aL = (r_a)e^iθ_1 , bL = (r_b)e^iθ_2 
        # Also find phi:
        # φ_L = θ_2 - θ_1

        # We must take into account the quadrants of aL and bL to find the correct polar angle θ, which is measured from the x-axis.
        '''
        S   |   A
        ____|____
            |   
        T   |   C
        
        '''

        if aLr > 0 and aLi >= 0:
            aL_theta = math.atan(abs(aLi)/abs(aLr))
        elif aLr < 0 and aLi >= 0:
            aL_theta = np.pi - math.atan(abs(aLi)/abs(aLr))
        elif aLr < 0 and aLi < 0:
            aL_theta = np.pi + math.atan(abs(aLi)/abs(aLr))
        elif aLr > 0 and aLi < 0:
            aL_theta = 2*np.pi - math.atan(abs(aLi)/abs(aLr))
        elif aLr == 0 and aLi > 0:
            aL_theta = np.pi/2
        elif aLr == 0 and aLi < 0:
            aL_theta = 3*np.pi/2
        elif aLr == 0 and aLi ==0:
            aL_theta = 0

        else:
            print("Didn't account for this αL = ",αL," in finding θ_L so have skipped this trajectory")
            print("aLi = ",aLi)
            print("aLr = ",aLr)
            aL_theta = 0
            filea = open("aL bugs, φ = {}π , θ = {}π.txt".format(φ/np.pi,θ/np.pi),"a")
            filea.write("Trajectory = ")
            filea.write(str(traj))
            filea.write("\n")
            filea.write("aLi = ")
            filea.write(str(aLi))
            filea.write("\n")
            filea.write("aLr = ")
            filea.write(str(aLr))
            filea.write("\n")
            filea.close()
            error = True

        bLr = βL.real
        bLi = βL.imag

        if bLr > 0 and bLi >= 0:
            bL_theta = math.atan(abs(bLi)/abs(bLr))
        elif bLr < 0 and bLi >= 0:
            bL_theta = np.pi - math.atan(abs(bLi)/abs(bLr))
        elif bLr < 0 and bLi < 0:
            bL_theta = np.pi + math.atan(abs(bLi)/abs(bLr))
        elif bLr > 0 and bLi < 0:
            bL_theta = 2*np.pi - math.atan(abs(bLi)/abs(bLr))
        elif bLr == 0 and bLi > 0:
            bL_theta = np.pi/2
        elif bLr == 0 and bLi < 0:
            bL_theta = 3*np.pi/2
        elif bLr == 0 and bLi ==0:
            bL_theta = 0
        else:
            print("Didn't account for this bL = ",βL,"so have skipped this trajectory")
            print("bLi = ",bLi)
            print("bLr = ",bLr)
            bL_theta = 0
            fileb = open("bL bugs, φ = {}π , θ = {}π.txt".format(phi/np.pi,th/np.pi),"a")
            fileb.write("Trajectory = ")
            fileb.write(str(traj))
            fileb.write("\n")
            fileb.write("bLi = ")
            fileb.write(str(bLi))
            fileb.write("\n")
            fileb.write("bLr = ")
            fileb.write(str(bLr))
            fileb.write("\n")
            fileb.close()
            error = True

    # print("φ = ",φ/np.pi,"π , θ = ",θ/np.pi,"π")

        # r_beta = math.sqrt((bLr**2) + (bLi**2))   #not used

        # use global phase to get a real only value for α (amplitude associated with state |0> on bloch sphere)
        phi = bL_theta - aL_theta #Want 0 ≤ ph < 2π
        # print("current φ = ",phi/np.pi,"π")

        #The current range that it's producing is 0 < th <= π and -2π <= phi <= 2π
        #Convert to domains for Rabi experiment:
            # -π/2 ≤ ph < π/2  (visually for great circles near XZ plane this makes more sense, rather than using 0 <= phi < π and having points with minutely negative y's being given a theta value of theta = 2π - theta on the Rabi curve)
            # 0 ≤ th < 2π (Rabi curve is |a_L|^2 vs. theta so want 0 ≤ th < 2π)
        pi = np.pi
        if phi == 2*np.pi or phi == -2*np.pi:
            phi = 0 
        if -2*np.pi < phi < -np.pi/2:
            phi = 2*np.pi + phi 
            #no need to change θ 
            #should now have -π/2 ≤ φ ≤ 2π 
        if 3*np.pi/2 <= phi < 2*np.pi:
            phi = phi - 2*np.pi
            #should now have -π/2 ≤ φ < 3π/2  
        if -np.pi/2 <= phi < np.pi/2:
            pass
        if np.pi/2 <= phi <= np.pi: #2nd quadrant inclusive
            phi = phi - np.pi 
            #this moves it from 2nd quadrant to -π/2 ≤ φ ≤ 0 
            theta = 2*np.pi - theta #reach the same point with an increased theta
        if np.pi < phi < 3*np.pi/2 :
            phi = phi - np.pi
            theta = 2*np.pi - theta
            #should now have -π/2 ≤ φ < π/2  
        if phi < -np.pi/2 or phi >= np.pi/2:
            print("φ out of expected domain, φ = ",phi/np.pi,"π")

        x = math.sin(theta)*math.cos(phi)
        y = math.sin(theta)*math.sin(phi)
        z = math.cos(theta)

    phiL = phi
    thetaL = theta  #(J'aurais dû l'appeler comme ça au début)

    return thetaL,phiL, x,y,z

def calcLogicalAlphaBeta(traj,contenta,contentb,a,b):

    aLeqn = contenta[traj] # This is the equation for alphaL (in terms of alpha (a) and beta (b) physical corresponding to the trajectory 
    bLeqn = contentb[traj] 

    # print(f"\naLeqn = {aLeqn}")
    # print(f"bLeqn = {bLeqn}")

    αLprime = eval(aLeqn) # This subs in a and b physical as calculated by giveMoreExactValues for this theta and phi
    βLprime = eval(bLeqn)

    threshold = 1e-8
    if abs(αLprime) < threshold:
        αLprime = 0

    if abs(βLprime) < threshold:
        βLprime = 0

    norm = np.sqrt(abs(αLprime)**2 + abs(βLprime)**2)

    flag = False
    
    if norm == 0:
        # print(f"      Unreachable trajectory {traj}")
        flag = True
        return 0, 0, flag
    
    # New code:
    elif αLprime == 0:
        return 0, 1, flag
    elif βLprime == 0:
        return 1, 0, flag
    # End of new code.


    else:

        # Normalise:
        αL = αLprime/norm
        βL = βLprime/norm
        




        # New code:
        if αL.imag < threshold:
            αL = αL.real
        
        # Make αL real and positive.

        if αL.real < 0:
            αL = - αL
            βL = - βL



        if αL.imag != 0:
            αLprime = αL * αL.conj()
            βLprime = βL * αL.conj()

            norm = np.sqrt(abs(αLprime)**2 + abs(βLprime)**2)

            αL = αLprime / norm
            βL = βLprime / norm

            if αL.imag < threshold:
                αL = αL.real
            if βL.imag < threshold:
                βL = βL.real

        # End of new code
        
        return αL, βL, flag

def rainbow_color_stops(n, end):
    return [ hls_to_rgb((end * i/(n-1)), 0.5, 1) for i in range(n) ]

def analytical_group_by_theta(thetaLarray, analyticalprobarray, xLp, yLp, zLp, colours):
    # This groups trajectories under the same theta (if they have the same theta)
    # and calculates the xLp, yLp, and zLp entries (cartesian coords) corresponding to each unique thetaL.
    # Also handles colours, which are expected to be RGB tuples.
    
    duplicates = defaultdict(list)
    # Create an enumerate object of the array (puts the index in front of each entry)
    enumerated_array = enumerate(thetaLarray)
    # Iterate over positions and numbers simultaneously (each count up in i also counts up in thetaL)
    for i, thetaL in enumerated_array:
        # Accumulate positions to the same number:
        duplicates[thetaL].append(i)  # duplicates is a dictionary with keys as thetaL and values as the positions in the array where that thetaL occurs.

    newthetas = list(duplicates.keys())
    positionslist = list(duplicates.values())

    num_thetas = len(newthetas)
    newanalyticalprobs = np.empty(num_thetas)
    newxLp = np.empty(num_thetas)
    newyLp = np.empty(num_thetas)
    newzLp = np.empty(num_thetas)

    # Initialize newcolours to store tuples, matching the data type of RGB tuples
    newcolours = np.empty(num_thetas, dtype=object)

    for i, theta_positions in enumerate(positionslist):
        # Use the first value at the index from theta_positions (assuming all values for the same theta are the same)
        newanalyticalprobs[i] = analyticalprobarray[theta_positions[0]]
        newxLp[i] = xLp[theta_positions[0]]
        newyLp[i] = yLp[theta_positions[0]]
        newzLp[i] = zLp[theta_positions[0]]
        newcolours[i] = colours[theta_positions[0]]  # Assign the RGB tuple directly

    newthetas = np.array(newthetas)  # Convert to numpy array
    
    return newthetas, newanalyticalprobs, newxLp, newyLp, newzLp, newcolours


# # # PLOTTER CODE:

if len(sys.argv) != 2:
    print("Input one argument (the code distance). E.g. python3 PLOT_ANALYTICAL_RABI.py 2")
    sys.exit(1)

distance = int(sys.argv[1])

rotation = 'rot'    # CHANGE ROTATED OR UNROTATED HERE

filea = open(f"logical_coefficients/aL_D{distance}_{rotation}ated.txt","r") 
fileb = open(f"logical_coefficients/bL_D{distance}_{rotation}ated.txt","r")
contenta = filea.readlines()
contentb = fileb.readlines()

# Replace ^ (put in by Mathematica) with ** (python exponentiation operator)
# Also remove { and } characters:
contenta = [line.replace('^', '**').replace('{', '').replace('}', '') for line in contenta]
contentb = [line.replace('^', '**').replace('{', '').replace('}', '') for line in contentb]

dist = distance
if rotation == 'unrot':
    nq = dist**2 + (dist-1)**2  # number of data qubits (unrotated surface code)
elif rotation == 'rot':
    nq = dist ** 2
ns = nq - 1                 # number of stabilisers
nt = 2**ns                  # number of trajectories

colours = rainbow_color_stops(nt, 0.9) # make as many colours as there are trajectories


thetasubdivisions = 256
phisubdivisions = 128


with open(f"Results_of_PLOT_ANALYTICAL_RABI/d{distance} {rotation} unreachable trajectories.txt", "a") as log_file, \
    open(f"Results_of_PLOT_ANALYTICAL_RABI/d{distance} {rotation} alpha_L.txt", "w") as alpha_file, \
    open(f"Results_of_PLOT_ANALYTICAL_RABI/d{distance} {rotation} beta_L.txt" , "w") as beta_file :

    ## Uncomment to do a range of phi's:
    # for i in range(phisubdivisions+1): 
        # physphicoef = i/phisubdivisions 
    
    ## Uncomment to just do a single phi:
    for i in [0]: 
        physphicoef = 0  # CHANGE INITIAL ROTATION OF ALL DATA QUBITS HERE !!


        ## Loop start:
        φ  = physphicoef * np.pi 


        ## Uncomment to do a single theta:
        for k in [0]: 
            physthetacoef = 1.3   # CHANGE INITIAL ROTATION OF ALL DATA QUBITS HERE !!
        
        ## Uncomment to do a range of thetas:
        # for k in range(thetasubdivisions+1):
            # physthetacoef = 2*(k/thetasubdivisions)
            

            ## Loop start:
            
            θ = physthetacoef*np.pi

            ## If the calculation's already been done, skip it:
            # file_path1 = "Results_of_PLOT_ANALYTICAL_RABI/rabi curves/d{}/θ={}π, φ={}π.png".format(distance, format(θ/np.pi, '.7f'), format(φ/np.pi, '.7f'))
            # file_path2 = "Results_of_PLOT_ANALYTICAL_RABI/bloch spheres/d{}/θ={}π, φ={}π.png".format(distance,format(θ/np.pi, '.7f'),format(φ/np.pi, '.7f'))
            # # Check if the file exists
            # if os.path.exists(file_path1):
            #     if os.path.exists(file_path2):
            #         print(f"d{distance} {rotation}, θ={physthetacoef:.7f}π, φ={physphicoef:.7f}π",end = '\r')
            #         continue

            print(f"d{distance} {rotation}, θ={physthetacoef:.7f}π, φ={physphicoef:.7f}π")#,end = '\r')
            
            alpha_file.write(f"d{distance} {rotation}, θ={physthetacoef:.7f}π, φ={physphicoef:.7f}π\n")#,end = '\r')

            beta_file.write(f"d{distance} {rotation}, θ={physthetacoef:.7f}π, φ={physphicoef:.7f}π\n")#,end = '\r')


            header_written = False # for the log file of unreachable trajectories

            a,b = giveMoreExactValues(θ,φ)
            print("a = cos(θ/2)      = ",a)
            print("b = sin(θ/2)*e^iφ = ",b)

            alpha_file.write(f"a = cos(θ/2)      = {a:.2f}\n")
            alpha_file.write(f"b = sin(θ/2)*e^iφ = {b:.2f} \n\n")
            
            beta_file.write(f"a = cos(θ/2)      = {a:.2f}\n")
            beta_file.write(f"b = sin(θ/2)*e^iφ = {b:.2f} \n\n")

            thetaLarray = np.empty(2**ns)
            thetaLarray[:] = np.nan
            analyticalprobarray = np.empty(2**ns)
            analyticalprobarray[:] = np.nan

            xLp = np.empty(2**ns)
            xLp[:] = np.nan
            yLp = np.empty(2**ns)
            yLp[:] = np.nan
            zLp = np.empty(2**ns)
            zLp[:] = np.nan

            for traj in range(nt):   
                αL, βL, flag = calcLogicalAlphaBeta(traj,contenta,contentb,a,b) # can get 'list index out of range' if manuall set a trajectory larger than the possible trajectories in the dist.


                # Printing for comparison with C code in git repo
                # print(f"{αL.real} {'+' if αL.imag >= 0 else '-'} {abs(αL.imag)} * I ,")
                # print(f"{βL.real} {'+' if βL.imag >= 0 else '-'} {abs(βL.imag)} * I ,")




                
                # Write traj and coefficients to coef_file:
                # alpha_file.write(bin(traj)[2:].zfill(ns) + "\n")

                alpha_file.write(f"{αL.real} {'+' if αL.imag >= 0 else '-'} {abs(αL.imag)} * I ,\n")
                beta_file.write(f"{βL.real} {'+' if βL.imag >= 0 else '-'} {abs(βL.imag)} * I ,\n")

                
                if flag == True:
                    # unreachable trajectory
                    if not header_written:
                        log_file.write(f"θ={physthetacoef}π, φ={physphicoef}π\n")
                        header_written = True
                    
                    # print("Casting unreachable trajectory", traj,"into the","\x1B[9mfire\x1B[0m" , "log file")

                    log_file.write(f"{traj},")
                    continue # go to next trajectory because this one's impossible
                
                analyticalprobarray[traj] = abs(αL)**2

                # Now calc. position on logical bloch sphere:
                thetaL,phiL,xL,yL,zL = calcLogicalCoords(αL,βL) 

                # print("θL = ",thetaL)
                # print("φL = ",phiL)

                if thetaL == xL == yL == zL == 0: continue # unreachable trajectory

                #Save analytical thetaL and points for trajectory which appeared in simulator. These are ordered by trajectory number, any blanks means the trajectory didn't appear in the simulation.

                thetaLarray[traj] = thetaL/np.pi   # (will be plotted in units of π)
                xLp[traj] = xL
                yLp[traj] = yL
                zLp[traj] = zL

            if header_written == True:
                log_file.write("\n")
            log_file.flush()


            # Above for loop has created:
            # - thetaLarray and analyticalprobarray, ordered by trajectory with analytical values for thetaL (in units of π) and abs(α)^2
            # - x,y,z points for visualising the analytical trajectories on the Bloch sphere


            # Each trajectory has a thetaL value. For trajectories with exactly the same thetaL, I will group their counts. groupedthetas is now an array of unique thetaL values and groupedprobs is the average of the probabilities of the grouped trajectories. 

            groupedthetas, groupedanalyticalprobs, groupedxLp, groupedyLp, groupedzLp, groupedcolours = analytical_group_by_theta(thetaLarray,analyticalprobarray,xLp,yLp,zLp,colours)

            # PLOTTING:
            titlefontsize = 25
            tickfontsize = 20
            labelfontsize = 25
            dotsize_analytical = 60  # (for scatter plot)
            dotsize_simulation = 60  #(for scatter plot)


            # Plot analytical Rabi curve:
            # add cos^2 curve:
            x = np.arange(0,2*np.pi,0.01*np.pi)   # start,stop,step
            y = (np.cos(x/2))**2 

            plt.figure(figsize=(9.6,7.2))
            plt.tight_layout()
            plt.plot(x/np.pi,y, c = 'C0',linewidth = 0.5)    #the cos^2 curve
            plt.xlabel('$θ_L$',fontsize = labelfontsize)
            plt.ylabel('P(0 logical)',fontsize = labelfontsize)
            plt.title('d{}{}, θ={}π, φ={}π'.format(distance,rotation,format(θ/np.pi, '.3f'),format(φ/np.pi, '.3f')),fontsize = titlefontsize)
            steps = [0,1/4,1/2,3*1/4,1,5*1/4,3*1/2,7*1/4,2*1]
            plt.xticks( steps, ['0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'], fontsize = tickfontsize)
            plt.yticks(fontsize = tickfontsize)

            # add analytical values of trajectories
            plt.scatter(groupedthetas,groupedanalyticalprobs, s = dotsize_analytical, color = groupedcolours,label = 'analytical')
            plt.legend(fontsize=tickfontsize)
            plt.tight_layout()

            plt.savefig("Results_of_PLOT_ANALYTICAL_RABI/rabi curves/d{}/θ={}π, φ={}π.png".format(distance,format(θ/np.pi, '.7f'),format(φ/np.pi,'.7f')))
            plt.savefig("Results_of_PLOT_ANALYTICAL_RABI/rabi curve.png")
            plt.clf()
            plt.close()

            # Bloch sphere analytical graph (visualising the analytical values on bloch sphere)

            B = qutip.Bloch()
            B.figsize = [5,6.2]
            B.make_sphere()
            B.zlabel = ["θ = {}π  , φ = {}π \n ∣0⟩ᴸ".format(format(round(θ/np.pi,2), '.2f'),format(round(φ/np.pi,2), '.2f')),'  ∣1⟩ᴸ'] 
            B.zlpos = [1.3,-1.25]
            pnts = [groupedxLp, groupedyLp, groupedzLp]
            B.point_size = [dotsize_analytical/1.5]
            B.point_marker = ["o"]
            B.point_color = ['#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) for r,g,b in groupedcolours]
            B.add_points(pnts,'m')
            B.render()
            B.fig.savefig("Results_of_PLOT_ANALYTICAL_RABI/bloch spheres/d{}/θ={}π, φ={}π.png".format(distance,format(θ/np.pi, '.7f'),format(φ/np.pi, '.7f')))
            B.fig.savefig("Results_of_PLOT_ANALYTICAL_RABI/bloch sphere.png")


            plt.clf()    #get rid of the analytical points (saved under curve - analytical)
            plt.close()

            filea.close()
            fileb.close()

cputime = (time.process_time() - debut)
executiontime = (time.time() - start)
print("")
print("")
print("CPU time = ",str(cputime),"s")
print("Real time = ",str(executiontime),"s")