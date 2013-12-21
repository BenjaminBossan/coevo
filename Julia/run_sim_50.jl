using DataFrames
using PyPlot

println("STARTED")


#println("random seed")
#srand(17413)


# P A R A M E T E R S

iterations=1::Int         #how often the simulation is repeated
gen=1::Int                  #number of simulated generations
w0=10.::Float64              #base fitness
tmax=50::Int                #time per generation
pA0=0.5::Float64            #initial success rate of option A
pB0=0.5::Float64            #initial success rate of option B
pskill=0.::Float64           #max var of (pos. or neg.) influence of skill
nindi=1000::Int            #total population
param=8::Int               #number of parameters per individual
q=.9::Float64               #discount rate for older memory
nstrat=10::Int              #number of strategies
kdoubt=2.::Float64           #threshold value for use of maj. tal. by IDCs
# environment
regime=1.::Float64           #how the environment changes. 0->no regression to the mean, 1->medium, 2->high
incr=0.02::Float64         #increment at which the environment becomes better or worse
pincr=1.0::Float64            #probability that environmental quality changes at all after each period
dpA=0.::Float64
dpB=0.::Float64              # changes in pA and pB;
pA0=pA0+dpA;pB0=pB0+dpB;


# INDICES
# for convenience, these indices are not passed to the coevo function.
# Changes here have thus to be updated in coevo as well
ichoice=1::Int;isucc=2::Int;ifit=3::Int;itally=4::Int;
ibias=5::Int;iskill=6::Int;istrat=7::Int;ibest=8::Int;

# F U N C T I O N S

# GENERATE ENVIRONMENT
function randomenvironment(tmax,regime,incr,pincr,pA0)

  # generate random vector for environmental change
  # determines whether there is a switch at all
  penvvecA=rand(1,tmax);
  penvvecA=penvvecA.<=pincr;
  # determines whether A or B increases
  envvecA=rand(1,tmax);

  # success rate of hunting grounds
  # initialization
  pA=zeros(1,tmax);
  pA[1]=pA0;

  # reversion parameter
  r=regime;

  # min and max p that can be attained
  minpA=pA0-((pA0-mod(pA0,incr))/incr)*incr;
  maxpA=pA0+((pA0-mod(pA0,incr))/incr)*incr;

  for t=2:tmax
      if penvvecA[t]==1
          if envvecA[t-1]<min(1,(max(0,.5-r*(pA[t-1]-.5))))
              pA[t]=min(pA[t-1]+incr,maxpA);
          else
              pA[t]=max(pA[t-1]-incr,minpA);
          end
      else
          pA[t]=pA[t-1];
      end
  end
  return round(pA, 2)
end

# LEARNING
function coevo50(
        tmax,    #tmax time of simulation
        nindi,   #num of individuals
        nstrat,  #num of strategy types
        incr,    #the increment
        pincr,   #probability that patch quality changes
        w0,      #the base fitness
        ninitial,#initial state of the population
        pA0,     #remember last pA
        pB0,     #remember last pB
        dpA,     #change in mean pA
        dpB,     #change in mean pB
        param,   #num of parameter values
        q,       #oblivousness, discount for older memory
        kdoubt,  #threshold value for ODCs
        regime);    #regime, how the environment changes. 1->low variance, 2->high variance

  # INDICES
  # for convenience, these indices are not passed to the coevo function.
  # Changes here have thus to be updated in coevo as well
  ichoice=1::Int;isucc=2::Int;ifit=3::Int;itally=4::Int;
  ibias=5::Int;iskill=6::Int;istrat=7::Int;ibest=8::Int;

  # INITIALIZE POPULATION TENSOR
  n=zeros(nindi,param,tmax+1);  #population tensor
  n[:,:,1]=ninitial;              #the initial state of the population
  # skill stays same
  n[:,iskill,:] = broadcast(+, ninitial[:,iskill], n[:,iskill,:])
  # strategy stays the same
  n[:,istrat,2:tmax] = broadcast(+, n[:,istrat,1], n[:, istrat, 2:tmax] )

  # change in mean pA and pB:
  # first subtract the change to get normal behavior of pA and pB
  pA0=pA0-dpA;pB0=pB0-dpB;
  pA = randomenvironment(tmax::Int,regime::Float64,
             incr::Float64,pincr::Float64,pA0::Float64)
  pB = randomenvironment(tmax::Int,regime::Float64,
             incr::Float64,pincr::Float64,pB0::Float64)

  # apply changes in mean pA and pB, consider the boundaries
  pA=min(1,max(0,pA+dpA));pB=min(1,max(pB+dpB,0));

  # vector containing the SKILL
  skillvec=ninitial[:,iskill];

  # random matrix used for determining success
  randsucmatA=rand(nindi,tmax);
  randsucmatB=rand(nindi,tmax);

  for t=1:tmax

    # random vector with size of population
    randIndi=rand(nindi,1);

    # S U C C E S S
    # update success
    n[:, isucc, t] = (
      (n[:, ichoice, t].==1).*(rand(nindi,1)-skillvec.<pA[t]) )
    n[:, isucc, t] += (
      (n[:, ichoice, t].==0).*(rand(nindi,1)-skillvec.<pB[t]) )

    # CHOICE FREQUENCIES
    # frequency of A choices in last period:
    x=sum(n[:,ichoice,t])/nindi;
    # freq of successful A choices
    xAs=sum((n[:,ichoice,t].==1)&(n[:,isucc,t].==1))/nindi;
    # freq of failed A choices
    xAf=sum((n[:,ichoice,t].==1)&(n[:,isucc,t].==0))/nindi;
    # freq of successful B choices
    xBs=sum((n[:,ichoice,t].==0)&(n[:,isucc,t].==1))/nindi;
    # freq of failed B choices
    xBf=sum((n[:,ichoice,t].==0)&(n[:,isucc,t].==0))/nindi;

    # update BIAS towards A or B for strategies relying on individual
    # learning, discount older bias
    n[:,ibias,t+1]=q.*n[:,ibias,t]+((n[:,ichoice,t].==n[:,isucc,t])*2-1)
    # bias+1 if A succeeded or B failed, -1 vice versa


    # F I T N E S S
    # the benefit vector
    # all successful individuals receive +1
    if t==1
      n[:,ifit,t]= w0 + n[:,isucc,t];
    else
      n[:,ifit,t] = (n[:,ifit,t-1]     # add fitness of last round
          + n[:,isucc,t] );            # add benefit if successful
    end

    # C H O I C E

    # the SCORE vector
    # this vector will determine which option an individual
    # chooses in the next round. A positive score leads to A
    # choice, a negative score to B choice.
    scvec=zeros(nindi,1);

    # 1 INDIDIVUAL LEARNING
    scvec += ( (n[:,istrat,t].==1).*           #if indiv. learning this round
        n[:,ibias,t+1] );                     #score corresponds to bias


    # 2 CONFORMISM

    scvec += (n[:,istrat,t].==2).*(-1 + 2*(randIndi.<((3-2*x)*x^2)));


    # 3 OPPORTUNISTIC INDIVIDUAL LEARNERS

    # if successful last period, use conformism
    index1 = ((n[:,istrat,t].==3)&(n[:,isucc,t].==1));
    scvec += index1.*(-1 + 2*(randIndi.<((3-2*x)*x^2)));

    # if unsuccessful last period, use individual learning
    scvec += ((n[:,istrat,t].==3)&(n[:,isucc,t].==0)).*n[:,ibias,t+1];


    # 4 OPPORTUNISTIC CONFORMISTS

    # if unsuccessful last period, use conformism
    index1 = ((n[:,istrat,t].==4)&(n[:,isucc,t].==0));
    scvec += index1.*(-1 + 2*(randIndi.<((3-2*x)*x^2)));

    # if successful last period, use individual learning
    scvec += ((n[:,istrat,t].==4)&(n[:,isucc,t].==1)).*n[:,ibias,t+1];


    # 5 IN DOUBT, CONFORM

    # if there is doubt, use conformism
    index1 = ((n[:,istrat,t].==5)&(abs(n[:,ibias,t]).<kdoubt));
    scvec += index1.*(-1 + 2*(randIndi.<((3-2*x)*x^2)));

    # if certain, use individual learning
    scvec += (n[:,istrat,t].==5).*(abs(n[:,ibias,t]).>=kdoubt).*n[:,ibias,t+1];


    # 6 IMITATE THE WEALTHIEST

    if sum(n[:,istrat,t].==6) > 0 #if there are ITW users at all

      #number of samples individuals
      talmax=3;

      #random sample of population
      randSample = ceil( nindi*rand(nindi, talmax) ); # indexes of random sample
      randSampleAddFit = rand(nindi, talmax)          # small random addition to fitness in order to break ties

      #choice , fitness of the random sample
      randSampleChoice = Int[]      # choice of random sample
      randSampleFit    = Int[]      # fitness of random sample
      randSampleSort   = Int[]      # sorted random sample
#      bestChoice = 0.
      vv = Float64[]
      ww = mean(n[:,ichoice,t])
      for ii = 1:nindi
        if n[ii, istrat] == 6
          randSampleFit    = n[ randSample[ii,:][:], ifit, t] + randSampleAddFit[ii, :][:]
          randSampleChoice = n[ randSample[ii,:][:], ichoice, t]
          randSampleSort   = sortperm(randSampleFit[:], rev=true)
          bestChoice       = randSampleChoice[randSampleSort[1]]
          scvec[ii]        = bestChoice
        end
      end

    end

    # 7 PBSL [4/-1]
    if sum(n[:,istrat,t].==7)>0
      scvec += ( (n[:,istrat,t].==7).*
          (-1+2*(randIndi.<xAs^7+21*xAf^5*xAs*(xAs+xBf)+7*xAs^6*(xBf+xBs)+21*xAs^5*(xBf+xBs)^2+35*xAs^4*(xBf+xBs)^3+xBf^6*(xBf+7*xBs)+21*xAs^2*xBf^3*(xBf^2+5*xBf*xBs+10*xBs^2)+(7*xAs*xBf^4*(2*xBf^2+12*xBf*xBs+15*xBs^2))/2+35*xAs^3*xBf*(xBf^3+4*xBf^2*xBs+6*xBf*xBs^2+4*xBs^3)+(35*xAf^4*xAs*(2*xAs^2+6*xBf^2+3*xAs*(2*xBf+xBs)))/2+35*xAf^3*(xAs^4+4*xAs*xBf^3+xBf^4+4*xAs^3*(xBf+xBs)+6*xAs^2*xBf*(xBf+2*xBs))+21*xAf^2*(xAs^5+xBf^5+5*xAs^4*(xBf+xBs)+10*xAs^3*(xBf+xBs)^2+10*xAs^2*xBf^2*(xBf+3*xBs)+5*xAs*xBf^3*(xBf+4*xBs))+7*xAf*(xAs^6+6*xAs^5*(xBf+xBs)+15*xAs^4*(xBf+xBs)^2+xBf^5*(xBf+3*xBs)+6*xAs*xBf^4*(xBf+5*xBs)+20*xAs^3*xBf*(xBf^2+3*xBf*xBs+3*xBs^2)+15*xAs^2*xBf^2*(xBf^2+4*xBf*xBs+6*xBs^2)))) );
    end

    # 8 PBSL [1/0]
    if sum(n[:,istrat,t].==8)>0
      scvec += ( (n[:,istrat,t].==8).*
          (-1+2*(randIndi.<(xAf^3+2*xAs^3+xBf^3+3*xAf^2*(2*xAs+xBf)+6*xAs^2*(xBf+xBs)+6*xAs*xBf*(xBf+xBs)+3*xAf*(2*xAs^2+xBf^2+2*xAs*(2*xBf+xBs)))/2)) );
    end

    # 9 PBSL McElreath
    if sum(n[:,istrat,t].==9)>0
      scvec += ( (n[:,istrat,t].==9).*
          (-1+2*(randIndi.<xAf^3+3*xAf^2*(xAs+xBf)+3*xAf*xAs*(xAs+2*xBf)+xAs*(xAs^2+3*xAs*(xBf+xBs)+3*xBf*(xBf+2*xBs)))) );
    end

    # 10 PBSL PAYOFF-CONFORMISM TRADE-OFF
    if sum(n[:,istrat,t].==10)>0
      scvec += ( (n[:,istrat,t].==10).*
          (-1+2*(randIndi.<xAf^6+6*xAf^5*(xAs+xBf)+15*xAf^4*(xAs+xBf)^2+10*xAf^3*(2*xAs^3+6*xAs^2*xBf+xBf^3+6*xAs*xBf*(xBf+xBs))+15*xAf^2*xAs*(xAs^3+4*xAs^2*(xBf+xBs)+6*xAs*xBf*(xBf+2*xBs)+2*xBf^2*(2*xBf+3*xBs))+6*xAf*xAs*(xAs^4+5*xAs^3*(xBf+xBs)+5*xBf^3*(xBf+2*xBs)+5*xAs^2*(2*xBf^2+4*xBf*xBs+xBs^2)+5*xAs*xBf*(2*xBf^2+6*xBf*xBs+3*xBs^2))+xAs*(xAs^5+6*xAs^4*(xBf+xBs)+15*xAs^3*(xBf+xBs)^2+6*xBf^3*(xBf^2+5*xBf*xBs+10*xBs^2)+10*xAs^2*(2*xBf^3+6*xBf^2*xBs+6*xBf*xBs^2+xBs^3)+15*xAs*xBf*(xBf^3+4*xBf^2*xBs+6*xBf*xBs^2+2*xBs^3)))) );
    end


    #score>0 -> choose A, else choose B
    n[:,ichoice,t+1] = max(sign(scvec),0);

#    println(size(n[:,ibest,t]), size(2*n[:,ichoice,t] -1) )
    n[:, ibest, t] = ( (2*n[:,ichoice,t] -1) .* sign(pA[t]-pB[t]) );
  end

  # cut last period
  n = n[:,:,1:tmax]

  return n, pA, pB

end

# NEXT GENERATION
function next_gen(nindi, param, nt, ichoice, ifit, itally, istrat)

  nNew=zeros(nindi,param);

  # extract relevant information
  ntf=nt[:,ifit];    #fitness
  ntc=nt[:,ichoice]; #choice
  ntt=nt[:,itally];  #sampling size
  nts=nt[:,istrat];  #strategy

  # REPLICATION
  # vector which will contain indices of offspring
  offspringind = Int[]

  # number of offspring proportional to fitness, queued
  for i=1:nindi
    append!(offspringind, i*ones(convert(Int, ntf[i])))
  end
  # shuffle offspring
  offspringind = offspringind[randperm(length(offspringind))]
  # next generation has same pop size
  offspring = offspringind[1:nindi]

  # populate the next generation
  nNew[:, ichoice] = ntc[offspring][:]
  nNew[:, itally]  = ntt[offspring][:]
  nNew[:, istrat]  = nts[offspring][:]

  return nNew
end

# I N I T I A L I Z A T I O N
ninitial=zeros(nindi,param);        #the initial state of the population
nt=zeros(nindi,param,gen);
if iterations<2
    perfmat=zeros(nstrat,gen);
elseif gen<2
    perfmat=zeros(nstrat,iterations);
end

# INITIAL POPULATION
ninitial[:,ichoice]=rand(nindi,1).>1/2;                #initialize random choice in 1st round
ninitial[:,iskill]=(-.5+rand(nindi,1))*pskill;        #skill is uniformly distributed with mean 0


# ~~~~~~~~~~~~ change initial population here~~~~~~~~~~~~~~~
# use this procedure to define the initial population for
# pure strategies
ind=[
    900;          # INDividual learners
    0;          # CONformists
    0;          # Opportunstic Individual Learners
    0;          # Opportunstic Conformists
    0;          # In Doubt, Conform
    100;          # Imitate The Wealthiest
    0;          # PBSLs [4/-1]
    0;          # PBSLs [1/0]
    0;          # PBSLs McElreath
    0];         # PBSLs Payoff-Conformism Trade-off

s = 0::Int
j = 1::Int
k = 0::Int
for numStrats in ind
  s+=1
  k+=numStrats
  ninitial[j:k, istrat] = s
  j+=numStrats
end

# unique strategies at the beginning
unikInit = unique(ninitial[:,istrat])

# measure duration of simulation
tic()
for x=1:iterations

  for g=1:gen

    # learning process
    (n, pA, pB) = coevo50(tmax,nindi,nstrat,incr,pincr,w0,
                  ninitial,pA0,pB0,dpA,dpB,param,q,kdoubt,regime)

    # LAST STATES
    choicesVec = n[:, ichoice, tmax] # remember last choices
    nt[:, :, g] = n[:, :, tmax] # remember final state
    pA0 = pA[1, tmax]
    pB0 = pB[1, tmax]

    # stats
    if iterations < 2
      #determine performance
      pident=find(pA==pB);
      nperf=n[:,ibest,:];
      nperf[:,pident]=[]; # remove draws
      nperf=[nperf+1]/2;
      unik = unique(n[:,istrat,1])
      for s = 1:length(unik)
        perfmat[unik[s], g] = mean(nperf[n[:,istrat,tmax].==unik[s],:])
      end
    end

#    ntn = n[:, :, tmax]

    ninitial2 = ninitial             # remember previous initial state
    # next generation
    nNew = next_gen(nindi, param, nt[:, :, g], ichoice, ifit, itally, istrat)
    ninitial[:, :] = nNew[: ,:]

    # check for extinction
    if mod(g, 100) == 0
      if length(unik) == 1
        println("All but strategy ", unik[1], " are extinct")
        break
      end
    end

    # progress
    if (mod(g, 50) == 0) & (gen>1) & (iterations==1)
      println("-- Progress: ", 100*g/gen, " %")
    end

  end

end
global n, pA, pB, unik

#Pkg.add("PyPlot")

#using("PyPlot")


# WRITE OUTPUT TO FILE

# one generation:
if (gen == 1) & (iterations == 1)
  randName = string("data/1gen" , string(ceil(10^6*rand()))[1:5] , ".csv")
  f = open(randName, "w")
  print(f, "pA;")
  print(f, "pB;")
  for s in unikInit
    print(f, "strat", convert(Int, s))
    if s != unikInit[length(unikInit)]
      print(f, ";")
    end
  end
  print(f, "\n")
  for t=1:tmax
    print(f, pA[t], ";")
    print(f, pB[t], ";")
    for s in unikInit
      print(f,
        sum(
          (n[:, istrat, t].==s)
          &
          (n[:, ichoice, t].==1)
        )/
        sum(
          (n[:, istrat, t].==s)
        )
      )
      if s != unikInit[length(unikInit)]
        print(f, ";")
      end
    end
    print(f, "\n")
  end
  close(f)
  # plot
  df1 = readtable(randName, separator=';');
  plot(df1[1]-df1[2]+.5, "k")
  for s = 3:length(unikInit)+2
    plot(df1[s])
  end
end

# several generations (evolution)
if (gen > 1) & (iterations == 1)
  randName = string("data/evo" , string(ceil(10^6*rand()))[1:5] , ".csv")
  f = open(randName, "w")
  # print header
  for s in unikInit
    print(f, "strat", convert(Int, s))
    if s != unikInit[length(unikInit)]
      print(f, ";")
    end
  end
  print(f, "\n")
  for gg = 1:gen
    for s in unikInit
      print(f,
        sum(
          (nt[:, istrat, gg].==s)
        )/nindi
      )
      if s != unikInit[length(unikInit)]
        print(f, ";")
      end
    end
    print(f, "\n")
  end
  close(f)
  # read dataframe and plot
  df1 = readtable(randName, separator=';');
  for s = 1:length(unikInit)
    plot(df1[s])
  end
end


toc()

println("FINISHED")


