def boundify(starName):
    ############################################################################
    #                           COLD BELT BOUNDS
    ############################################################################
    ###############
    # UPPER BOUNDS
    ###############
    if starName == 'HD 110411':
        # Indicative of cleaner ice water in the outer region since it needs
        # a narrower SED for the cold belt
        blsBoundUpperOut = 1
    elif starName == 'HD 113337':
        blsBoundUpperOut = 5

    else: # default is 5
        blsBoundUpperOut = 5

    # Lower bound is how much the bound is divided by. Use a fraction
    # to increase the lower bound
    ###############
    # LOWER BOUNDS
    ###############
    if starName == 'HD 110897':
        blsBoundLowerOut = 1/3
    elif starName == 'HD 113337':
        blsBoundLowerOut = 1
    else: # default is 10
        blsBoundLowerOut = 10


    ############################################################################
    #                           WARM BELT BOUNDS
    ############################################################################
    ###############
    # UPPER BOUNDS
    ###############
    if starName == 'HD 110411':
        blsBoundUpperIn = 3
    else:
        blsBoundUpperIn = 3

    ###############
    # LOWER BOUNDS
    ###############
    # Lower bound is how much the bound is divided by. Use a fraction
    # to increase the lower bound
    if starName == 'HD 110411':
        blsBoundLowerIn = 10
    else:
        blsBoundLowerIn = 10


    return blsBoundLowerIn, blsBoundUpperIn, blsBoundLowerOut, blsBoundUpperOut
