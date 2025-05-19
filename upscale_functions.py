def gpp_model(Age, model, T=None, Vpd=None, Fert=None, Habitat=0):
    """fixed part of gpp-age model
    Args:
        Age (float): stand age (yr)
        model (int): model
        T (float): mean air temperature (degC). Defaults to None.
        Vpd (float): growing season vapor pressure deficit (kPa). Defaults to None.
        Fertility (int): Low == 1, High = 0. Defaults to None.
        Habitat (int): Hemiboreal == 1, Boreal = 0
 
    Returns:
        float: modeled GPP (g C m-2 a-1)
    """
   
    # model 1-4 parameterized with all data!
    if model == 1: # Age-only
        a = -71.67
        b = 1493.49
        c = 6.92
   
    if model == 2: # Age, T
        a = 14.60
        b = 861.21 + 92.79 * T
        c = 5.45
 
    if model == 3:  # Age, T, VPD_gs
        a = 52.210274    
        b = 1223.00 + 86.70 * T - 414.48 * Vpd
        c = 5.54
 
    if model == 4: # Age, T, VPD_gs, Fertility
        a = 52.210274    
        b = 1223.00 + 86.70 * T - 414.48 * Vpd - 422.61 * Fert
        c = 5.54
 
    # model 5-8 parameterized with subset of data, neglecting DK, IS sites, SE-Nor and part of EE sites!
    if model == 5: # Age-only
        a = -84.23
        b = 1427.76
        c = 6.07
   
    if model == 6: # Age, T
        a = -12.20
        b = 1000.22 + 62.84 * T
        c = 5.22
 
    if model == 7:  # Age, T, VPD_gs
        a = 0.51    
        b = 1084.58 + 67.52 * T - 292.95 * Vpd
        c = 5.28
 
    if model == 8: # Age, T, VPD_gs, Fertility
        a = 34.07    
        b = 1228.71 + 42.66 * T - 246.60 * Vpd - 403.72 * Fert + 337.32 * Habitat
        c = 5.54
 
    y = a + b * (Age / (Age + c))
 
    return y
 
def reco_model(Age, model, T=None, Plant=None, Fert=None):
    """fixed part of reco-age model
    Args:
        Age (float): stand age (yr)
        model (int, optional): model
        T (float): mean air temperature (degC). Defaults to None.
        Plant (int): Plantation = 1, other = 0
        Fertility (int): Low == 1, High = 0. Defaults to None.
       
    Returns:
        float: modeled RECO (g C m-2 a-1)
    """    
    # Full dataset used in parameterization    
    if model == 1: # Age-only
        a = 878.95
        d = 1.69
       
    if model == 2: # Age, T
        a = 687.70 + 38.97 * T
        d = 1.53
       
    if model == 3:
        a = 950.76 + 31.02 * T - 496.52 * Fert - 327.00 * Plant    
        d = 2.14
    # parameterized with subset of data, neglecting DK, IS sites, SE-Nor and part of EE sites!
   
    if model == 4: # Age-only
        a = 925.25
        d = 0.77
       
    if model == 5: # Age, T
        a = 789.14 + 26.98 * T
        d = 0.78
       
    if model == 6:
        a = 993.87 + 25.63 * T - 512.26 * Fert  
        d = 1.86
 
    y = a + d * Age
 
    return y
 
def nep_model(Age, model, T=None, Plant=None, Fert=None):
    """fixed part of nep-age model
    Args:
        Age (float): stand age (yr)
        model (int, optional): model. Defaults to 4.
        T (float): mean air temperature (degC). Defaults to None.
        Plant (int): Plantation = 1, other = 0
        Fert (int): Low == 1, High = 0. Defaults to None.
       
    Returns:
        float: modeled NEP (g C m-2 a-1)
    """  
    # Full dataset used in parameterization    
    if model == 1: # Age-only
        a = -639.11
        b = 1100.45
        c = 7.17
        d = 1.70
 
    if model == 2: # Age, T
        a = -606.20
        b = 897.62 +  31.85 * T
        c = 7.14
        d = 1.84
   
    if model == 3:
        a =  -632.19 + 405.46 * Plant
        b = 843.81 + 32.88 * T
        c = 6.56
        d = 1.54
   
    # parameterized with subset of data, neglecting DK, IS sites, SE-Nor and part of EE sites!
 
    if model == 4: # Age-only
        a = -638.52
        b = 1086.75
        c = 6.81
        d = 1.90
 
    if model == 5: # Age, T
        a = -613.73
        b = 891.00 +  31.64 * T
        c = 6.73
        d = 1.91
   
    if model == 6: # Age, T, Fertility
        a =  -632.61 + 58.69 * Fert
        b = 924.47 + 35.27 * T
        c = 6.45
        d = 1.87
 
    if model == 6: # Age, T, Fertility, Habitat
        a =  -697.46 + 99.57 * Fert + 141.88 * Habitat
        b = 852.55 + 26.76 * T
        c = 6.95
        d = 2.23
 
    y = a + b * (Age / (Age + c)) - d * Age
 
    return y