


# packages
import copy
import numpy  as np
import pandas as pd
import scipy
import scipy.interpolate



# =============================================================================
# Implementation of ModVege [ Jouven et al. (2006) ] - DM changes directly linked to LAI values | water and nitrogen liminted production
# =============================================================================


### model inputs
# (i)   - proportion of functional grass groups;
# (ii)  - N nutrion index;
# (iii) - soil water holding capacity;
# (iv)  - weather data [ incoming radiation, temperature (max and min), precipitation and potential evapotranspiration ];
# (v)   - cutting regime;
# (vi)  - sward status in the begining of the simulation period ; approximate values can be given (results for 1st year are discarted) or
#         model is run twice (1st to estimate initial biomass/sward age at a given date and 2nd to produce final estimates) - to be checked/improved



#### test
#ClimData     = ClimDataSumm
#HarvestDates = HarvestDatesL
#NIval        = NIvalL[ NIvalL.index(0.85) ]
#FieldLabel   = FieldLabel
#SoilPropL    = SoilPropL
#WaterResFuncOpt = 'Option1'
#EnvFuncOpt      = 'Option1'
#GroundWContr    = 0.2





###
class ModVege:


# =============================================================================
# initialization
# =============================================================================

 #
 def __init__( self, swardComp, ClimData ):

  ## parameters common to all functional groups (to be changed if necessary)
  self.sigmGveg = 0.4     # rate of biomass loss with respiration - Green vegetative component
  self.sigmGrep = 0.2     # rate of biomass loss with respiration - Green reprodutive component
  self.TT0      = 3.0     # base temp. - 4.0 deg (original) , LINGRA = 3.0 deg 
  self.TT1      = 10.0    # optimal temp.
  self.TT2      = 20.0    # limit of optimal temp. - stress starting
  self.KGveg    = 0.002   # senescence rate - Green vegetative
  self.KGrep    = 0.001   # senescence rate - Green reprodutive
  self.KlDveg   = 0.001   # abscission rate - Dead vegetative
  self.KlDrep   = 0.0005  # abscission rate - Dead reprodutive
  self.OMDDveg  = 0.45    # organic matter digestibility - Dead vegetative material
  self.OMDDrep  = 0.40    # organic matter digestibility - Dead reprodutive material

  ## parameters according to sward composition - Cruz et al. (2002) ; (values for each functional groups - to be changed if necessary)
  self.sla        = np.multiply( [ 0.033, 0.025, 0.022, 0.019 ], swardComp ).sum().round(4)  # specific leaf area [ m2 g-1 ]
  self.lam        = np.multiply( [ 0.680, 0.680, 0.680, 0.680 ], swardComp ).sum().round(4)  # percentage of laminae in the accumulated material
  self.st1        = np.multiply( [ 600, 700, 850, 1000 ]       , swardComp ).sum().round(4)  # thermal time = initial growth
  self.st2        = np.multiply( [ 1200, 1350, 1550, 1850 ]    , swardComp ).sum().round(4)  # thermal time = end of growth
  self.SEAmax     = np.multiply( [ 1.20, 1.30, 1.40, 1.50 ]    , swardComp ).sum().round(4)  # max seasonal effect
  self.SEAmin     = np.multiply( [ 0.80, 0.70, 0.60, 0.50 ]    , swardComp ).sum().round(4)  # min seasonal effect
  self.lls        = np.multiply( [ 500, 800, 900, 1400 ]       , swardComp ).sum().round(4)  # leaf life span
  self.OMDmaxGveg = np.multiply( [ 0.90, 0.90, 0.85, 0.75 ]    , swardComp ).sum().round(4)  # max organic material digestibility (OMD) - Green vegetative
  self.OMDminGveg = np.multiply( [ 0.75, 0.60, 0.65, 0.65 ]    , swardComp ).sum().round(4)  # min OMD - Green vegetative
  self.OMDmaxGrep = np.multiply( [ 0.90, 0.90, 0.85, 0.75 ]    , swardComp ).sum().round(4)  # max OMD - Green reprodutive
  self.OMDminGrep = np.multiply( [ 0.65, 0.45, 0.45, 0.45 ]    , swardComp ).sum().round(4)  # min OMD - Green reprodutive
  self.BDGveg     = np.multiply( [ 850, 850, 1200, 800 ]       , swardComp ).sum().round(4)  # bulk density (BD - g DM m-3) - Green vegetative
  self.BDDveg     = np.multiply( [ 500, 500, 1800, 2200 ]      , swardComp ).sum().round(4)  # BD - Dead vegetative
  self.BDGrep     = np.multiply( [ 300, 300, 200, 150 ]        , swardComp ).sum().round(4)  # BD - Green reproductive
  self.BDDrep     = np.multiply( [ 150, 150, 300, 450 ]        , swardComp ).sum().round(4)  # BD - Dead reproductive

  ## other parameters used in the model (not related to the functional groups ; to be changed if necessary)
  self.k         = 0.6   # canopy light extinction coefficient
  self.RUEmax    = 3.0   # maximum radiation use efficiency [ g DM MJ-1 ]; same value used in LINGRA-N ; 2.2-2.4 g DM / MJ of PAR = Insua et al., 2019
  self.kc        = 1.15  # crop coef. for actual ET calculation
  self.PARfrac   = 0.47  # fraction of the incoming radiation corresponfing to the PAR (Photosynthetically Active Radiation)
  self.multFact  = 15    # multiplicative factor used to extend interval of valid inputs for functions describing crop development
  self.reduFact  = 0.01  # multiplicative factor used to create small deviation between values in fuctions with non-linear response


#  # minimum levels allowed for biomass in different components
#  self.NullBMGvegThrsld = 0
#  self.NullBMGrepThrsld = 0
#  self.NullBMDvegThrsld = 0
#  self.NullBMDrepThrsld = 0


  ## model outputs - creation of data frame in which daily steps will be stored
  self.ModDayOut = pd.DataFrame( np.nan, index=ClimData.index, columns=['BMGveg', 'BMGrep', 'BMDveg', 'BMDrep', \
                                                                        'AGEGveg', 'AGEGrep', 'AGEDveg', 'AGEDrep', \
                                                                        'OMDGveg', 'OMDGrep', 'LAI', 'waterContRel', \
                                                                        'EvapTranspCrop', 'PAR', 'SEAeff', 'fAPAR', 'CumTemp', \
                                                                        'REPpart', 'REPindi', 'stressWater', 'stressTemp', \
                                                                        'stressPAR', 'EnvFactEnsamb' ] )



# =============================================================================
# Defining model functions
# =============================================================================


 ## 1) Plant response to daily temperature - as implemented by Schapedonk et al., 1998
 def plantRespT(self, DayAvgTemp ):
  Resp_outT = np.max( [0, np.min( [1, ((DayAvgTemp - self.TT0) / (self.TT1 - self.TT0))] )] )
  return Resp_outT


 ## 2) Plant response to PAR level
 def plantRespPAR(self, DayPAR ):
  Resp_outT = np.min( [1, 1-0.0445*(DayPAR-5)] )
  if Resp_outT < 0:
   Resp_outT = 0
  return Resp_outT


 ## seasonal effects
 ## 3) Sum of temperature since 1st January of the year considered
 def TempCumYear(self, DateStart, DateCurrent, ClimData, FieldTempMax, FieldTempMin ):
  ## DateStart = starting date ; DateCurrent = current date in Timestamp format ; ClimData = climatic data as pandas dataframe with index in Timestamp format
  ## FieldTempMax, FieldTempMin = str() indicating the variable names containing max and min temperature in the climatic dataframe
  # calculate average temperature for each day in the period
  TempMaxPer = ClimData.loc[ DateStart:DateCurrent, FieldTempMax ].values
  TempMinPer = ClimData.loc[ DateStart:DateCurrent, FieldTempMin ].values
  TempAvgPer = ( (TempMaxPer+TempMinPer)/2 ).tolist()
  TempAvgPerO = [ TempAvgPer[val]-self.TT0 if TempAvgPer[val] > self.TT0 else 0 for val in range(len(TempAvgPer)) ]   # only consider temperature above base
  TempAvgPerO = np.nansum( TempAvgPerO )  # get thermal time
  if TempAvgPerO < 0:
   TempAvgPerO = 0
  return TempAvgPerO

 ## 4) Function describing seasonal effects according to TempCumYear
 def SEAeff(self, DateStart, DateCurrent, ClimData, FieldTempMax, FieldTempMin ):
  # get cumulative temperature since 1st January of the same year
  TempCumYearO = self.TempCumYear( DateStart, DateCurrent, ClimData, FieldTempMax, FieldTempMin )
  # values used to describe seasonal effect (linear interpolation between values)
  SEAval = [ self.SEAmin, self.SEAmin,  self.SEAmax,  self.SEAmax, self.SEAmin, self.SEAmin ]
  TTval  = [           0,         200, self.st1-200, self.st1-100,    self.st2, self.multFact*self.st2 ]
  SEAeff_funcT = scipy.interpolate.interp1d( TTval, SEAval,  kind='linear' )
  return SEAeff_funcT( TempCumYearO ).tolist()


 ## 5) Reproductive function
 def REPf(self, NIval, DateStart, DateCurrent, ClimData, FieldTempMax, FieldTempMin ):
  ## NIval = Nutrition Index [ between 0.35 and 1.20 according to Calanca et al. (2016) ; indicate nitrogen status ]
  TempCumYearO = self.TempCumYear( DateStart, DateCurrent, ClimData, FieldTempMax, FieldTempMin )
  CUTO    = 1 if all([ (TempCumYearO >= self.st1), (TempCumYearO <= self.st2) ]) else 0  # boolean function indicating if reprodutive growth occurs
  REPoutT = ( 0.25 + ( (( 1 - 0.25 ) * ( NIval - 0.35 )) / ( 1 - 0.35 ) ) ) * CUTO
  return REPoutT


 ## 6) Plant resp to age - Senescence - green vegative
 def plantRespAgeSENGveg(self, AGEGveg ):
  ## AGEGveg = average age of green vegetative compartiment in thernal time (same in functions for other compartiments)
  AGEGvegT                 = [ 0, (1/3)*self.lls, (2/3)*self.lls, self.lls, self.multFact*self.lls ]
  plantRespAgeSENGveg_outT = [ 1, 1, 2, 3, 3 ]
  plantRespAgeSENGvegT     = scipy.interpolate.interp1d( AGEGvegT,
                                                         plantRespAgeSENGveg_outT, kind='linear' )
  return plantRespAgeSENGvegT( AGEGveg ).tolist()


 ## 7) Plant resp to age - Senescence - green reproductive
 def plantRespAgeSENGrep(self, AGEGrep ):
  AGEGrepT                 = [ 0, (1/3)*(self.st2-self.st1), (2/3)*(self.st2-self.st1), (self.st2-self.st1), self.multFact*(self.st2-self.st1) ]
  plantRespAgeSENGrep_outT = [ 1, 1, 2, 3, 3 ]
  plantRespAgeSENGrepT     = scipy.interpolate.interp1d( AGEGrepT,
                                                         plantRespAgeSENGrep_outT, kind='linear' )
  return plantRespAgeSENGrepT( AGEGrep ).tolist()


 ## 8) Plant response to age - Abscission - dead vegetative
 def plantRespAgeABSDveg(self, AGEDveg ):
  AGEDvegT                 = [ 0, ((1/3)*self.lls)-self.reduFact, ((1/3)*self.lls),
                                  ((2/3)*self.lls)-self.reduFact, ((2/3)*self.lls),
                               self.lls, self.multFact*self.lls ]
  plantRespAgeABSDveg_outT = [ 1, 1, 2, 2, 3, 3, 3 ]
  plantRespAgeABSDvegT = scipy.interpolate.interp1d( AGEDvegT,
                                                     plantRespAgeABSDveg_outT, kind='linear' )
  return plantRespAgeABSDvegT( AGEDveg ).tolist()


 ## 9) Plant response to age - Abscission - dead reproductive
 def plantRespAgeABSDrep(self, AGEDrep ):
  AGEDrepT                 = [ 0, ((1/3)*(self.st2-self.st1))-self.reduFact, ((1/3)*(self.st2-self.st1)),
                                  ((2/3)*(self.st2-self.st1))-self.reduFact, ((2/3)*(self.st2-self.st1)),
                               (self.st2-self.st1), self.multFact*(self.st2-self.st1) ]
  plantRespAgeABSDrep_outT = [ 1, 1, 2, 2, 3, 3, 3 ]
  plantRespAgeABSDrep = scipy.interpolate.interp1d( AGEDrepT,
                                                    plantRespAgeABSDrep_outT, kind='linear' )
  return plantRespAgeABSDrep( AGEDrep ).tolist()


 ## 10) Water holding capacity [ m3 m-3 ] derived from soil properties (clay, sand, organic matter and root zone depth)
 def sWHC(self, clayC, sandC, OrMatterC, depthRZ ):  # depthRZ = depth of root zone [ meters ]
  OrMatterCc = OrMatterC * 0.001 * 10  # transform [ g kg-1 ] >> [ kg kg 1:1000 ] >> [ % ]
  # -33 kPa - field capacity
  fieldCapT = 0.2576 + (-0.002*sandC) + (0.0036*clayC) + (0.0299*OrMatterCc)
  # -1500  kPa - wilting point
  wiltPT = 0.026 + (0.005*clayC) + (0.0158*OrMatterCc)
  #
  deltaCapT = fieldCapT - wiltPT
  holdCapT  = np.round( [deltaCapT * depthRZ * 1000], 4 )[0]  # [ mm m-3 ] of soil  ;  x1000 to get per 1 m-3 ; x800 to get per 0.8 m-3  [ Lazzaroto et al., 2009 ]
  return holdCapT


 ## 11) Plat response to water | opt1 - FAO 56 document [ Allen et al., 1998 ] - smoother
 def plantRespW_opt1(self, wT, potEvapTransp ):  # wT = ratio between actual and potential water content / soil holding capacity
  p5T = 0.6  # for grasslands
  pT = p5T + 0.04*(5-potEvapTransp)
  wCrT = 1-pT
  OutT = np.min( [(wT/wCrT),1] )
  return OutT


 ## 12) Plat response to water | opt2 - Calanca et al. (2016) ; McCall and Bishop-Hurley (2003)
 def plantRespW_opt2(self, wT, potEvapTransp ):  # wT = ratio between actual and potential water content / soil holding capacity
  #
  if potEvapTransp <= 3.81:
   if wT <= 0.2:
    OutT = 4*wT
   else:
    if all([(wT > 0.2), (wT <= 0.4)]):
     OutT = 0.75*wT + 0.65
    else:
     if all([(wT > 0.4), (wT <= 0.6)]):
      OutT = 0.25*wT + 0.85
     else:
      OutT = 1
  #
  else:
   if all([(potEvapTransp > 3.81), (potEvapTransp <= 6.35)]):
    if wT <= 0.2:
     OutT = 2*wT
    else:
     if all([(wT > 0.2), (wT <= 0.4)]):
      OutT = 1.5*wT + 0.1
     else:
      if all([(wT > 0.4), (wT <= 0.6)]):
       OutT = wT + 0.3
      else:
       if all([(wT > 0.6), (wT <= 0.8)]):
        OutT = 0.5*wT + 0.6
       else:
        OutT = 1
   #
   else:
    OutT = wT
  return OutT



# =========================================================
# run growth model - iteration over all years and cut dates
# =========================================================


 ## this function runs the model according to the input data
 def runModel( self, ClimDataIn, HarvestDates, NIval, FieldLabel, SoilProp, WaterResFuncOpt, EnvFuncOpt, GroundWContr, minDurRepr, verboseOpt ):
  ## inputs for this function should be:
  # ClimDataIn      = pandas dataframe with climatic data;
  # HarvestDates    = list with harvest dates in Timestamp format;
  # NIval           = the nitrogen index level to consider;
  # FieldLabel      = list with the labels for each respective field in the climatic dataset [ 'MaxTemp', 'MinTemp', 'ETo', 'InRadiation', 'Precip' ];
  # SoilProp        = list with the soil properties for the site of interest to calculate the soil holding capacity [ clay, sand, organic matter and root zone depth ];
  # WaterResFuncOpt = function used to describe the plant response to water content [ 'Option1' or 'Option2' ]
  # EnvFuncOpt      = function used to describe aggregate environmental factors constraint [ 'Option1' or 'Option2' ]
  # GroundWContr    = minimum water content possible in the soil profile [ 0 to 1 w.r.t soil holding capacity ]
  # minDurRepr      = minimum duration (days) of the reproductive stage

  #
  YearRefAll = np.unique( [HarvestDates[val].year for val in range(len(HarvestDates))] ).tolist()
  for ii in range( len(YearRefAll) ):  # year to be simulated
   # calculations starting for 1st January of each year ; one interval added to cover each complete year (i.e., to always reach 31st December)
   DateStartSel = pd.Timestamp( int( YearRefAll[ii] ), int( 1 ), int( 1 ) )  # starting date for the year (1st January)

   #
   HarvYear = np.array(HarvestDates)[ [(HarvestDates[val].year == YearRefAll[ii]) for val in range(len(HarvestDates))] ]  # dates for cuts in the seleted year
   LastCut  = pd.Timestamp( int( YearRefAll[ii] ), int( 12 ), int( 31 ) )  # add last period for the given year - simulate the complte year from 1-Jan tp 31-Dec (remember to not harvest this date...)
   HarvYear = np.append( HarvYear, np.array([LastCut]), axis=0 )

   # check if time-series is complete for a given year
   if sum(ClimDataIn.index.year == YearRefAll[ii]) not in [365, 366]:
    print( 'Climatic data for ' + str(YearRefAll[ii]) + ' is not complete (or DataFrame is not correctly indexed)' )
    break


   #
   for jj in range( HarvYear.shape[0] ):  # harvest dates in the selected year
    if ( verboseOpt == True ):
     print( 'Simulating harvest for:', HarvYear[jj].strftime("%d-%m-%Y"), sep=' ' )


    # select current period ( from previous cut to current harvest date )
    if (jj==0):  # 1st harvest of the year
     PrevCut = pd.Timestamp( int( YearRefAll[ii] ), int( 1 ), int( 1 ) )
    else:
     PrevCut = HarvYear[jj-1]
    # Dates in the period to be processed
    PeriodDates = ClimDataIn.loc[ PrevCut:HarvYear[jj],: ].index.tolist()


    #
    for kk in range( len(PeriodDates) ):  # daily steps for the period between current and previous cut
     DateCurrentSel = PeriodDates[kk]     # current date


     ### model initialization
     if all([ (ii==0), (jj==0), (kk==0) ]):  # 1st year, harvest and day
      # initial water content ; considering that initial water content is maximal | can be replaced by previously estimated value in the future (?)
      waterCont = 1.0  # [ % / 100 ]

      # residual material left in the field after previous harvest | assuming that residual is always the same (5 cm hight stalks) | can be improved in the future (?)
      resBMGveg = 0.05 * 10 * self.BDGveg  # 420
      resBMGrep = 0                        # 0
      resBMDveg = 0.05 * 10 * self.BDDveg  # 300
      resBMDrep = 0                        # 30
      # age for each residual component
      AGEGveg = 400 # 100
      AGEGrep = 0   # 2000
      AGEDveg = 500 # 500
      AGEDrep = 0   # 500
      #
      BMGveg = copy.deepcopy( resBMGveg )  # not included in the outputs (only after one day of growth...)
      BMGrep = copy.deepcopy( resBMGrep )
      BMDveg = copy.deepcopy( resBMDveg )
      BMDrep = copy.deepcopy( resBMDrep )
      # these values will be used to assign current BM levels | OMD in the second 'IF' statement below (1st day after each cut), if this is the first year simulated

      # digestibility - only to check in this case...
      OMDGveg = self.OMDmaxGveg - ( ( AGEGveg * ( self.OMDmaxGveg - self.OMDminGveg ) ) / self.lls )
      OMDGrep = self.OMDmaxGrep - ( ( AGEGrep * ( self.OMDmaxGrep - self.OMDminGrep ) ) / ( self.st2 - self.st1 ) )

      #
      LAIcur = self.sla * (BMGveg/10) * self.lam


     #
     if all([ (jj==0), (kk==0) ]): # 1st harvest and day for each year
      # indicators for reproductive function - one cycle per year (0 = not started yet)
      repOccurCur  = 0
      repOccurNext = 0


     #
     if all([ (jj!=0), (kk==0) ]):  # 1st day after each cut ( except after givenYear-12-31 , when harvest does not occur - previous values propagated through new iteration )

      ## initialization assuming the same proportion for each compartiment in the residual biomass | vegetation 5cm long (after harvest)
      # green vegetative
      if all([ (BMGveg >= (0.05 * 10 * self.BDGveg)) ]):  # current biomass content longer than 5 cm (assuming a given bulk density)
       resBMGveg = 0.05 * 10 * self.BDGveg
      else:  # biomass shorter than 5 cm (not removed)
       if all([ (BMGveg > 0), (BMGveg < (0.05 * 10 * self.BDGveg)) ]):
        resBMGveg = copy.deepcopy(BMGveg)
       else:  # no biomass in the compartiment [ zero value propagated to next cycle - however in this case growth would not occur (no leaves left to start accumulation...) ]
        resBMGveg = 0
      # dead vegetative
      if all([ (BMDveg >= (0.05 * 10 * self.BDDveg)) ]):
       resBMDveg = 0.05 * 10 * self.BDDveg
      else:
       if all([ (BMDveg > 0), (BMDveg < (0.05 * 10 * self.BDDveg)) ]):
        resBMDveg = copy.deepcopy(BMDveg)
       else:
        resBMDveg = 0
      # green reproductive
      if all([ (BMGrep >= (0.05 * 10 * self.BDGrep)) ]):
       resBMGrep = 0.05 * 10 * self.BDGrep
      else:
       if all([ (BMGrep > 0), (BMGrep < (0.05 * 10 * self.BDGrep)) ]):
        resBMGrep = copy.deepcopy(BMGrep)
       else:
        resBMGrep = 0
      # dead reproductive
      if all([ (BMDrep >= (0.05 * 10 * self.BDDrep)) ]):
       resBMDrep = 0.05 * 10 * self.BDDrep
      else:
       if all([ (BMDrep > 0), (BMDrep < (0.05 * 10 * self.BDDrep)) ]):
        resBMDrep = copy.deepcopy(BMDrep)
       else:
        resBMDrep = 0

      #
      BMGveg = copy.deepcopy( resBMGveg )  # not included in the outputs (only after one day of growth)
      BMGrep = copy.deepcopy( resBMGrep )
      BMDveg = copy.deepcopy( resBMDveg )
      BMDrep = copy.deepcopy( resBMDrep )

      #
      LAIcur = self.sla * (BMGveg/10) * self.lam

      # propagate indicator of the reproductive function to the next period
      repOccurCur = copy.deepcopy( repOccurNext )


     # save bulk density values for the current sward composition (calculated as a function of sward composition)
     self.ModDayOut.loc[ PrevCut:HarvYear[jj], 'BDGveg' ] = copy.deepcopy( self.BDGveg )
     self.ModDayOut.loc[ PrevCut:HarvYear[jj], 'BDGrep' ] = copy.deepcopy( self.BDGrep )
     self.ModDayOut.loc[ PrevCut:HarvYear[jj], 'BDDveg' ] = copy.deepcopy( self.BDDveg )
     self.ModDayOut.loc[ PrevCut:HarvYear[jj], 'BDDrep' ] = copy.deepcopy( self.BDDrep )


     ### average air temperature | daily calculation
     tempCur = ( (ClimDataIn.loc[DateCurrentSel,FieldLabel[0]] + ClimDataIn.loc[DateCurrentSel,FieldLabel[1]]) / 2 )
#     if tempCur < self.self.TT0:  # if temperature is below base (?)
#      tempCur = 0


     ### calc water stress coefficient for the current day - proportional to ( water content / soil holding capacity ) - Calanca et al., 2016
     sWHCr = self.sWHC( SoilProp[0], SoilProp[1], SoilProp[2], SoilProp[3] )
     wT = (( waterCont * sWHCr ) / sWHCr)  # available water vs field capacity - not used directly as water stress factor since relationship is not linear ( plantRespW used instead ) ; in this case multiplication and division by sWHCr is not necessary, used just to indicate that value is given w.r.t. sWHCr
     # water stress factor | plant response to wT status - daily
     if (WaterResFuncOpt == 'Option1'):
      feW = self.plantRespW_opt1( wT,  (ClimDataIn.loc[DateCurrentSel,FieldLabel[2]] * self.kc) )
     else:
      feW = self.plantRespW_opt2( wT,  (ClimDataIn.loc[DateCurrentSel,FieldLabel[2]] * self.kc) )
     # save values - daily
     self.ModDayOut.loc[ DateCurrentSel, 'stressWater' ] = copy.deepcopy( feW )


     # potential  Transp | Evap
     potTransp = ( ClimDataIn.loc[DateCurrentSel,FieldLabel[2]] * self.kc ) * ( 1-np.exp(-self.k*LAIcur) )
     potEvap   = ( ClimDataIn.loc[DateCurrentSel,FieldLabel[2]] * self.kc ) - potTransp
     # actual  Transp | Evap
     actTransp = potTransp * feW  # correct plant transpiration according to plant response to water status (depending on Pot ET)
     actEvap   = potEvap * wT
     ETcCur    = actTransp + actEvap  # current ETc - daily
     # save values - daily
     self.ModDayOut.loc[ DateCurrentSel, 'EvapTranspCrop' ] = copy.deepcopy( ETcCur )


     ### incoming radiation
     PARcur = ClimDataIn.loc[DateCurrentSel,FieldLabel[3]]*self.PARfrac
     # save values - daily
     self.ModDayOut.loc[ DateCurrentSel, 'PAR' ] = copy.deepcopy( PARcur )


     ### seasonal effect for the day
     SEAcur = self.SEAeff( DateStartSel, DateCurrentSel, ClimDataIn, FieldLabel[0], FieldLabel[1] )
     # save values - daily
     self.ModDayOut.loc[ DateCurrentSel, 'SEAeff' ] = copy.deepcopy( SEAcur )


     ### plant response to air temperature | according to Schapedonk et al., 1998  or  Jouven et al., 2006 (change response function)
     feT = self.plantRespT( tempCur )
     # save values - daily
     self.ModDayOut.loc[ DateCurrentSel, 'stressTemp' ] = copy.deepcopy( feT )


     ### plant response to PAR
     fePAR = self.plantRespPAR( PARcur )
     # save values - daily
     self.ModDayOut.loc[ DateCurrentSel, 'stressPAR' ] = copy.deepcopy( fePAR )


     ### emsamble of environmental factors
     if (EnvFuncOpt == 'Option1'):
      # opt 1 - original
      ENVens = NIval * feW * fePAR * feT
     else:
      # opt 2 - consider nitrogen and water as competing stress fractors (min value retained) - according to [ Bonesmo et al., 2002 ]
      ENVens = np.nanmin( [NIval, feW] ) * fePAR * feT  # consider nitrogen and water as competing stress fractors (min value retained)
     # save values - daily
     self.ModDayOut.loc[ DateCurrentSel, 'EnvFactEnsamb' ] = copy.deepcopy( ENVens )


     ### fraction of Absorbed PAR
     fAPARcur = 1 - np.exp( (-self.k*LAIcur) )
     # save values - daily
     self.ModDayOut.loc[ DateCurrentSel, 'fAPAR' ] = copy.deepcopy( fAPARcur )


     ### potential growth
     PGROcur = PARcur * self.RUEmax * fAPARcur * 10  # [ kg ha-1 ] - daily


     ### actual growth
     GROcur = PGROcur * ENVens * SEAcur  # [ kg ha-1 ] - daily


     ### reproductive growth - one cycle per year
     if (repOccurCur == 0):
      REPcur = self.REPf( NIval, DateStartSel, DateCurrentSel, ClimDataIn, FieldLabel[0], FieldLabel[1] )
      # indicate that the reproductive phase was activated in the current cut ; at least 10 days of reproductive growth (check in the literature a better criteria)
      if all([ (repOccurCur == 0), \
               (self.TempCumYear( DateStartSel, DateCurrentSel, ClimDataIn, FieldLabel[0], FieldLabel[1] ) >= self.st1), \
               ((PeriodDates[-1] - DateCurrentSel).days >= minDurRepr) ]):
       repOccurNext = 1
     else:
      REPcur = 0
     # save values  - reproductvive partitioning - daily
     self.ModDayOut.loc[ DateCurrentSel, 'CumTemp' ] = copy.deepcopy( self.TempCumYear( DateStartSel, DateCurrentSel, ClimDataIn, \
                                                                                        FieldLabel[0], FieldLabel[1] ) )
     self.ModDayOut.loc[ DateCurrentSel, 'REPpart' ] = copy.deepcopy( REPcur )
     #
     self.ModDayOut.loc[ DateCurrentSel, 'REPindi' ] = copy.deepcopy( repOccurCur )



     ### green vegetative growth - daily
     GROGveg = GROcur * (1-REPcur)

     ### green reproductive growth - daily
     GROGrep = GROcur * REPcur # daily


     ### Senescence - vegetative
     if ( tempCur > self.TT0 ):
      SENGveg = self.KGveg * BMGveg * (tempCur) * self.plantRespAgeSENGveg( AGEGveg )  # as described in the methodlogy ( not using (tempCur - self.TT0) )
     else:
      if ( tempCur < 0 ):  # freezing effect
       SENGveg = self.KGveg * BMGveg * abs(tempCur)
      else:
       SENGveg = 0


     ### Senescence - reproductive | one cycle per year
     if ( tempCur > self.TT0 ):
      SENGrep = self.KGrep * BMGrep * (tempCur) * self.plantRespAgeSENGrep( AGEGrep )  # as described in the methodlogy ( not using (tempCur - self.TT0) )
     else:
      if ( tempCur < 0 ):  # freezing effect
       SENGrep= self.KGrep * BMGrep * abs(tempCur)
      else:
       SENGrep = 0


     ### rate of biomass change - green vegetative material
     rBMGveg = GROGveg - SENGveg

     ### rate of biomass change - green reproductive material
     rBMGrep = GROGrep - SENGrep


     ### abscission
     if ( tempCur > 0 ):
      ABSDveg = self.KlDveg * BMDveg * (tempCur) * self.plantRespAgeABSDveg( AGEDveg )
      ABSDrep = self.KlDrep * BMDrep * (tempCur) * self.plantRespAgeABSDrep( AGEDrep )
     #
     else: # abscission occuring only if temperature is above 0dC
      ABSDveg = 0
      ABSDrep = 0


     ### rate of biomass change (dead vegetative or reproductive material)
     rBMDveg = (1-self.sigmGveg) * SENGveg - ABSDveg
     rBMDrep = (1-self.sigmGrep) * SENGrep - ABSDrep


     ### age for each component
     if ( tempCur > 0 ):
      if  (( BMGveg - SENGveg + GROGveg ) != 0):
       rAGEGveg = ( ( BMGveg - SENGveg ) / ( BMGveg - SENGveg + GROGveg ) ) * ( AGEGveg + tempCur ) - AGEGveg
      else:
       rAGEGveg = 0
      #
      if  (( BMGrep - SENGrep + GROGrep ) != 0):
       rAGEGrep = ( ( BMGrep - SENGrep ) / ( BMGrep - SENGrep + GROGrep ) ) * ( AGEGrep + tempCur ) - AGEGrep
      else:
       rAGEGrep = 0
      #
      if (( BMDveg - ABSDveg + (1-self.sigmGveg) * SENGveg ) != 0):
       rAGEDveg  = ( ( BMDveg - ABSDveg ) / ( BMDveg - ABSDveg + (1-self.sigmGveg) * SENGveg ) ) * ( AGEDveg + tempCur ) - AGEDveg
      else:
       rAGEDveg = 0
      #
      if (( BMDrep - ABSDrep + (1-self.sigmGrep) * SENGrep ) != 0):
       rAGEDrep = ( ( BMDrep - ABSDrep ) / ( BMDrep - ABSDrep + (1-self.sigmGrep) * SENGrep ) ) * ( AGEDrep + tempCur ) - AGEDrep
      else:
       rAGEDrep = 0


     ### update pools | param.
     # biomass
     BMGveg = BMGveg + rBMGveg
#     if BMGveg < self.NullBMGvegThrsld:
#      BMGveg = self.NullBMGvegThrsld
     self.ModDayOut.loc[ DateCurrentSel, 'BMGveg' ] = copy.deepcopy( BMGveg ) # save values - daily
     #
     BMGrep = BMGrep + rBMGrep
#     if BMGrep < self.NullBMGrepThrsld:
#      BMGrep = self.NullBMGrepThrsld
     self.ModDayOut.loc[ DateCurrentSel, 'BMGrep' ] = copy.deepcopy( BMGrep ) # save values - daily
     #
     BMDveg = BMDveg + rBMDveg
#     if BMDveg < self.NullBMDvegThrsld:
#      BMDveg = self.NullBMDvegThrsld
     self.ModDayOut.loc[ DateCurrentSel, 'BMDveg' ] = copy.deepcopy( BMDveg ) # save values - daily
     #
     BMDrep = BMDrep + rBMDrep
#     if BMDrep < self.NullBMDrepThrsld:
#      BMDrep = self.NullBMDrepThrsld
     self.ModDayOut.loc[ DateCurrentSel, 'BMDrep' ] = copy.deepcopy( BMDrep ) # save values - daily


     # age
     AGEGveg = AGEGveg + rAGEGveg
     self.ModDayOut.loc[ DateCurrentSel, 'AGEGveg' ] = copy.deepcopy( AGEGveg ) # save values - daily
     #
     AGEGrep = AGEGrep + rAGEGrep
     self.ModDayOut.loc[ DateCurrentSel, 'AGEGrep' ] = copy.deepcopy( AGEGrep ) # save values - daily
     #
     AGEDveg = AGEDveg + rAGEDveg
     self.ModDayOut.loc[ DateCurrentSel, 'AGEDveg' ] = copy.deepcopy( AGEDveg ) # save values - daily
     #
     AGEDrep = AGEDrep + rAGEDrep
     self.ModDayOut.loc[ DateCurrentSel, 'AGEDrep' ] = copy.deepcopy( AGEDrep ) # save values - daily


     # digestibility - vegetative [ constant for dead components ]
     OMDGveg = self.OMDmaxGveg - ( ( AGEGveg * ( self.OMDmaxGveg - self.OMDminGveg ) ) / self.lls )
     self.ModDayOut.loc[ DateCurrentSel, 'OMDGveg' ] = copy.deepcopy( OMDGveg ) # save values - daily
     #
     OMDGrep = self.OMDmaxGrep - ( ( AGEGrep * ( self.OMDmaxGrep - self.OMDminGrep ) ) / ( self.st2 - self.st1 ) )
     self.ModDayOut.loc[ DateCurrentSel, 'OMDGrep' ] = copy.deepcopy( OMDGrep ) # save values - daily


     # update water status | LAI after daily growth
     precSel    = ClimDataIn.loc[DateCurrentSel,FieldLabel[4]]
     DeltaInOut = precSel-ETcCur
     waterCont  = ( (waterCont + (DeltaInOut / sWHCr)) if ( (waterCont + (DeltaInOut / sWHCr)) < (sWHCr/sWHCr) ) else (sWHCr/sWHCr) )  # waterCont given w.r.t. sWHCr
     if (waterCont < 0):
      waterCont = 0
     #
     if (GroundWContr > 0):  # indicate if ground water contribution is assumed to occur - limiting minimum water content observed in the profile
      waterCont = copy.deepcopy(GroundWContr) if all([ (waterCont < GroundWContr), (GroundWContr <= 1) ]) else waterCont # replace current water content by the minimum threshold allowed
     #
     self.ModDayOut.loc[ DateCurrentSel, 'waterContRel' ] = copy.deepcopy( waterCont ) # save values - daily


     # update LAI according to biomass accumulated
     LAIcur   = self.sla * (BMGveg/10) * self.lam
     self.ModDayOut.loc[ DateCurrentSel, 'LAI' ] = copy.deepcopy( LAIcur ) # save values - daily


