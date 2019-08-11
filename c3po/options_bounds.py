'''
------------Set the bounds for individual stars using this file-----------------

################################################################################
     Template for changing the bound options. These are the defaults!
################################################################################


    They must always be in ascending order, from Low -> High.
    HighInner/HighOuter must always be a greater number than LowInner/LowOuter


elif starName == 'HD 105':
    # Warm
    amin1_Low  = 0.1
    amin1_High = 3
    # Cold
    amin2_Low  = 0.1
    amin2_High = 3



'''
def boundify(starName):
    if starName == 'HD 105':
        # Warm
        amin1_Low  = 0.1
        amin1_High = 2
        # Cold
        amin2_Low  = 3
        amin2_High = 5
    elif starName == 'HD 10647':
        # Warm
        amin1_Low  = 0.1
        amin1_High = 2
        # Cold
        amin2_Low  = 1.25
        amin2_High = 5
    elif starName == 'HD 32297':
        # Warm
        amin1_Low  = 0.1
        amin1_High = 0.7
        # Cold
        amin2_Low  = 0.1
        amin2_High = 2
    elif starName == 'HD 110897':
        # Warm
        amin1_Low  = 0.1
        amin1_High = 3
        # Cold
        amin2_Low  = 0.1
        amin2_High = 3
    elif starName == 'HD 15115':
        # Warm
        amin1_Low  = 0.1
        amin1_High = 3
        # Cold
        amin2_Low  = 1
        amin2_High = 3

    elif starName == 'HD 153053':
        # Warm
        amin1_Low  = 0.1
        amin1_High = 3
        # Cold
        amin2_Low  = 1
        amin2_High = 3
    elif starName == 'HD 164249':
        # Warm
        amin1_Low  = 0.1
        amin1_High = 3
        # Cold
        amin2_Low  = 1
        amin2_High = 3
    elif starName == 'HD 192425':
        # Warm
        amin1_Low  = 0.1
        amin1_High = 3
        # Cold
        amin2_Low  = 0.01
        amin2_High = 0.5
    # These are the default values
    else:
        # Warm
        amin1_Low  = 0.1
        amin1_High = 3
        # Cold
        amin2_Low  = 0.1
        amin2_High = 3

    return amin1_Low, amin1_High, amin2_Low, amin2_High
