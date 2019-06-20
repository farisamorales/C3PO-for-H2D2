'''
Set the bounds for individual stars using this file.
'''

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
        bosHighOuter = 1
    elif starName == 'HD 113337':
        bosHighOuter = 5

    else: # default is 5
        bosHighOuter = 5

    # Lower bound is how much the bound is divided by. Use a fraction
    # to increase the lower bound
    ###############
    # LOWER BOUNDS
    ###############
    if starName == 'HD 110897':
        bosLowOuter = 1/3
    elif starName == 'HD 113337':
        bosLowOuter = 1
    else: # default is 10
        bosLowOuter = 10


    ############################################################################
    #                           WARM BELT BOUNDS
    ############################################################################
    ###############
    # UPPER BOUNDS
    ###############
    if starName == 'HD 110411':
        bosHighInner = 3
    else:
        bosHighInner = 3

    ###############
    # LOWER BOUNDS
    ###############
    # Lower bound is how much the bound is divided by. Use a fraction
    # to increase the lower bound
    if starName == 'HD 110411':
        bosLowInner = 10
    else:
        bosLowInner = 10


    return bosLowInner, bosHighInner, bosLowOuter, bosHighOuter
