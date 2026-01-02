"""
Truss in uss - Team 806
"""

"""
center.txt
10 0 12 73.92
88 0 90 73.92
10 73.92 15 74.92
85 73.92 90 74.92
0 74.92 100 80

sides.txt
10 0 12 75.19
88 0 90 75.19
10 75.19 15 76.19
85 75.19 90 76.19
0 76.19 100 80

ends.txt
10 0 12 75.19
88 0 90 75.19
10 75.19 15 76.19
85 75.19 90 76.19
12 0 88 1.27
0 76.19 100 80

design_zero.txt
10 0 11.27 73.73
11.27 0 88.73 1.27
10 73.73 16.27 75
88.73 0 90 73.73
83.73 73.73 90 75
0 75 100 76.27
"""

import matplotlib.pyplot as plt
import math

def sfd(load_setup, front_pos):
    # Determine positions of loads on the beam
    load_pos = []
    for i in range(6):
        pos = front_pos - train_load_loc[i]
        if 0 <= pos <= 1200:
            load_pos.append((pos, load_setup[i]))
    
    load_pos.sort(key=lambda x: x[0])
    
    # Calculating reactions
    rb = 0
    sum_loads = 0
    for pos, load in load_pos:
        if 0 <= pos <= 1200:
            rb += load * pos / 1200
            sum_loads += load
    ra = sum_loads - rb
        
    # Constructing SFD
    shear = ra
    sfd = [(0, shear)]
    for pos, load in load_pos:
        if sfd[-1][0] < pos:
            sfd.append((pos, shear))
        
        shear -= load
        
        if sfd[-1][0] == pos:
            sfd[-1] = (pos, shear)
        else:
            sfd.append((pos, shear))
    
    # Ending the SFD at 0
    if sfd[-1][0] != 1200:
        sfd.append((1200, 0))
    else:
        sfd[-1] = (1200, 0)

    return sfd

def sfd_envelope(load_setup, step_size):
    max = [-float('inf')]*1201
    min = [float('inf')]*1201
    sfd_matrix = []
    
    # Simulates moving the train across the bridge in increments of step_size
    for front_loc in range(0, 1200 + 176*3 + 164*2 + 1, step_size):
        sfd_i = sfd(load_setup, front_loc)
                
        sfd_data = []
        x = 0
        sfd_idx = 0
        current_shear_for_plot = sfd_i[0][1]
        
        # Filling in the SFD data at every position along the span
        while sfd_idx < len(sfd_i):
            while x < sfd_i[sfd_idx][0]:
                sfd_data.append((x, current_shear_for_plot))
                x += 1
            current_shear_for_plot = sfd_i[sfd_idx][1]
            sfd_data.append((sfd_i[sfd_idx][0], current_shear_for_plot))
            sfd_idx += 1
        
        # Updating the max and min envelopes
        for pos, shear in sfd_data:
            if shear > max[int(pos)]:
                max[int(pos)] = shear
            if shear < min[int(pos)]:
                min[int(pos)] = shear
        
        sfd_matrix.append(sfd_data)
        
        # Plotting each SFD
        plt.plot([pos for pos, shear in sfd_data], [shear for pos, shear in sfd_data])
    
    # Plot setup
    plt.title("All Shear Force Diagrams")
    plt.xlabel("Position Along Span (mm)")
    plt.ylabel("Shear Force (N)")
    plt.grid(True, linestyle="--")
    plt.plot(range(1201), max, color='red', linewidth=2, label='Envelope (max)')
    plt.plot(range(1201), min, color='blue', linewidth=2, label='Envelope (min)')
    plt.xlim(0,1200)
    # plt.show()
    
    return sfd_matrix, max, min

def bmd_envelope(load_setup, step_size):
    max = [-float('inf')]*1201
    
    # Using SFD data to do Riemann sum for BMD
    sfd_matrix, max_sfd, min_sfd = sfd_envelope(load_setup, step_size)
    
    for sfd in sfd_matrix:
        bmd_i = []
        moment = 0
        for pos, shear in sfd:
            moment += shear
            bmd_i.append((pos, moment))
        
        for pos, moment in bmd_i:
            if max[pos] < moment:
                max[pos] = moment
        
        # Plotting each BMD
        plt.plot([pos for pos, moment in bmd_i], [moment for pos, moment in bmd_i])
    
    # Plot setup
    plt.title("All Bending Moment Diagrams")
    plt.xlabel("Position Along Span (mm)")
    plt.ylabel("Moment (Nmm)")
    plt.grid(True, linestyle="--")
    plt.plot(range(1201), max, color='red', linewidth=2, label='Envelope')
    plt.xlim(0,1200)
    ax = plt.gca()
    ax.invert_yaxis()

    # Display the plot
    # plt.show()
    
    return max

def cross_sectional_properties():
    global centroidal_axis, centroid_width, cross_section_height
    global I
    global Q
    global Q_glue, glue_tab_width
    
    # location of the centroidal axis from the bottom
    rectangle_coords = []
    with open(filename, 'r') as f:
        for line in f:
            floats = [float(x) for x in line.split()]
            left_bottom = (floats[0], floats[1])
            right_upper = (floats[2], floats[3])
            rectangle_coords.append((left_bottom, right_upper))

    numer = 0
    denom = 0
    cross_section_height = 0
    y = 0
    for i in range(len(rectangle_coords)):
        base = rectangle_coords[i][1][0] - rectangle_coords[i][0][0]
        height = rectangle_coords[i][1][1] - rectangle_coords[i][0][1]
        if rectangle_coords[i][1][1] > y:
            y = rectangle_coords[i][1][1]
            cross_section_height += height
        y_i = (height/2) + rectangle_coords[i][0][1]
        A_i = base * height
        numer += A_i*y_i
        denom += A_i

    centroidal_axis = numer/denom
    print("Centroidal Axis: " + str(centroidal_axis) + " mm")


    # I (Second Moment of Area)
    I = 0
    for i in range(len(rectangle_coords)):
        base = rectangle_coords[i][1][0] - rectangle_coords[i][0][0]
        height = rectangle_coords[i][1][1] - rectangle_coords[i][0][1]
        y_i = (height/2) + rectangle_coords[i][0][1]
        A_i = base * height
        d_i = abs(centroidal_axis - y_i)
        first_term = (base * height**3)/12
        second_term = A_i*d_i**2
        I += first_term + second_term
    print("Second Moment of Area: " + str(I) + " mm^4")


    # Q at centroidal axes
    # we can first look at every single piece to see if it's over the centroid
    i = 0
    rectangle_coords_copy = rectangle_coords.copy()
    while i < len(rectangle_coords_copy):
        if rectangle_coords_copy[i][0][1] < centroidal_axis:
            # now let's check to see if it's fully below or if it intersects
            if rectangle_coords_copy[i][1][1] > centroidal_axis:
                # accumulating the total centroid width
                centroid_width += abs(rectangle_coords_copy[i][1][0] - rectangle_coords_copy[i][0][0])
                
                # then we can change it to end at the centroid
                temp = ((rectangle_coords_copy[i][0][0], rectangle_coords_copy[i][0][1]), (rectangle_coords_copy[i][1][0], centroidal_axis))
                rectangle_coords_copy[i] = temp
            i += 1
        else:
            rectangle_coords_copy.pop(i)
    Q = 0
    for i in range(len(rectangle_coords_copy)):
        base = rectangle_coords_copy[i][1][0] - rectangle_coords_copy[i][0][0]
        height = rectangle_coords_copy[i][1][1] - rectangle_coords_copy[i][0][1]
        y_i = (height/2) + rectangle_coords_copy[i][0][1]
        A_i = base * height
        d_i = centroidal_axis - y_i
        Q += A_i*d_i
    print("First moment of area at centroid: " + str(Q) + " mm^3")
    
    # Q at glue location (assume only the deck is above the glue spots)
    glue_height = float(input("Please input bottom y coordinate of glue tabs: "))
    glue_tab_width = rectangle_coords[2][1][0] - rectangle_coords[2][0][0]
    base = rectangle_coords[-1][1][0] - rectangle_coords[-1][0][0]
    height = rectangle_coords[-1][1][1] - rectangle_coords[-1][0][1]
    A_i = base * height
    d_i = abs(rectangle_coords[-1][0][1] + height/2 - centroidal_axis)
    Q_glue = A_i*d_i
    print("First Moment of Area at Glue: " + str(Q_glue) + " mm^3")

def max_shear(load_setup, step_size):
    sfd_matrix, maximum, min = sfd_envelope(load_setup, step_size)
    
    # shear is highest at the ends of the span of each cross-section, so it just returns the higher of the two
    if filename == 'sides.txt':
        print("Max Shear: " + str(max(maximum[70], abs(min[70]), maximum[1200-70], abs(min[1200-70]))) + " Nmm")
        return max(maximum[70], abs(min[70]), maximum[1200-70], abs(min[1200-70]))
    elif filename == 'center.txt':
        print("Max Shear: " + str(max(maximum[400], abs(min[400]), maximum[1200-400], abs(min[1200-400]))) + " Nmm")
        return max(maximum[400], abs(min[400]), maximum[1200-400], abs(min[1200-400]))
    else:
        print("Max Shear: " + str(max(round(maximum[0]), round(abs(min[-2])))) + " Nmm")
        return max(round(maximum[0]), round(abs(min[-2]))) 

def max_moment(load_setup, step_size):
    max_bmd = bmd_envelope(load_setup, step_size)
    
    # checking for close-to-maximum bending moments for different cross sections
    if filename == 'sides.txt':
        print("Max Moment: " + str(max(max_bmd[300], max_bmd[1200-300])) + " Nmm")
        return max(max_bmd[300], max_bmd[1200-300])
    elif filename == 'ends.txt':
        print("Max Moment: " + str(max(max_bmd[50], max_bmd[1200-50])) + " Nmm")
        return max(max_bmd[50], max_bmd[1200-50])
    max_moment_value = 0
    for moment in max_bmd:
        if moment > max_moment_value:
            max_moment_value = moment
    
    print("Max Moment: " + str(max_moment_value) + " Nmm")
    return max_moment_value

def applied_stresses(load_setup, step_size):
    max_moment_value = max_moment(load_setup, step_size)
    max_shear_value = max_shear(load_setup, step_size)
    flex_top = (max_moment_value * (cross_section_height - centroidal_axis)) / I
    flex_bottom = (max_moment_value * (centroidal_axis)) / I
    shear_centroid = (max_shear_value * Q) / (I * centroid_width)
    shear_glue = (max_shear_value * Q_glue) / (I * glue_tab_width * 2)

    return flex_top, flex_bottom, shear_centroid, shear_glue

def thin_plate_buckling():
    global centroidal_axis
    global I
    mu = 0.2
    E = 4000
    # getting coordinates from cross-section file
    rectangle_coords = []
    with open(filename, 'r') as f:
        for line in f:
            floats = [float(x) for x in line.split()]
            left_bottom = (floats[0], floats[1])
            right_upper = (floats[2], floats[3])
            rectangle_coords.append((left_bottom, right_upper))
            
    # let's find vertical and horizontal members (ASSUME ONLY TWO VERTICAL)
    vertical = []
    horizontal = []
    for i in range(len(rectangle_coords)):
        # we're looking for vertical members because their height > base
        base = rectangle_coords[i][1][0] - rectangle_coords[i][0][0]
        height = rectangle_coords[i][1][1] - rectangle_coords[i][0][1]
        if base > height:
            horizontal.append(rectangle_coords[i])
        else:
            vertical.append(rectangle_coords[i])

    # now let's look for top deck members
    max_height = 0
    deck_length = 0
    for i in range(len(horizontal)):
        if horizontal[i][1][1] > max_height:
            max_height = horizontal[i][1][1]
            deck_length = horizontal[i][1][0] - horizontal[i][0][0]
    deck_pieces = []
    for i in range(len(horizontal)):
        length = horizontal[i][1][0] - horizontal[i][0][0]
        if length == deck_length:
            deck_pieces.append(horizontal[i])

    # let's find the thickness of the top deck
    deck_thickness = 0
    for i in range(len(deck_pieces)):
        deck_thickness += (deck_pieces[i][1][1] - deck_pieces[i][0][1])

    # let's find the thickness of one vertical member
    vertical_thickness = vertical[1][1][0] - vertical[1][0][0]

    # let's find the thickness of the bottom member
    bottom_thickness = 0
    for i in range(len(horizontal)):
        if horizontal[i][0][1] == 0:
            bottom_thickness = horizontal[i][1][1] - horizontal[i][0][1]

    # Case 1:
    t1 = deck_thickness # thickness of the top deck
    if filename == 'design_zero.txt':
        b1 = horizontal[2][0][0] - horizontal[1][1][0]
    else:
        b1 = deck_length - (horizontal[1][0][0] - horizontal[0][1][0])
    S1 = (((4*(math.pi)**2)*E)/(12*(1-mu**2)))*(t1/b1)**2

    # Case 2:
    t2 = deck_thickness # thickness of the top deck
    b2 = min(vertical[1][0][0], vertical[0][0][0])# distance from the end to the first vertical member
    S2 = ((0.425*((math.pi)**2)*E)/(12*(1-mu**2)))*(t2/b2)**2

    # Case 3:
    t3 = vertical_thickness # thickness of the vertical
    b3 = max_height - deck_thickness - centroidal_axis # distance between global centroidal and bottom of the deck
    S3 = ((6*((math.pi)**2)*E)/(12*(1-mu**2)))*(t3/b3)**2


    # Case 4:
    t4 = vertical_thickness # thickness of the vertical pieces
    h4 = max_height - bottom_thickness - deck_thickness # the height until the top flange
    a4 = diaphragm_gap # distance between the center of each vertical member
    Tau4 = ((5*((math.pi)**2)*E)/(12*(1-mu**2)))*((t4/h4)**2 + (t4/a4)**2)

    print(S1, S2, S3, Tau4)
    return S1, S2, S3, Tau4

def fos(stresses, thin_plate_buckling):
    flex_top, flex_bottom, shear_centroid, shear_glue = stresses
    S1, S2, S3, Tau4 = thin_plate_buckling
    flex_compression_allowable = 6 
    flex_tension_allowable = 30
    shear_allowable = 4
    shear_glue_allowable = 2

    fos_tension = flex_tension_allowable / flex_bottom
    fos_compression = flex_compression_allowable / flex_top
    fos_shear_centroid = shear_allowable / shear_centroid
    fos_shear_glue = shear_glue_allowable / shear_glue
    fos_buckling_1 = S1 / flex_top
    fos_buckling_2 = S2 / flex_top
    fos_buckling_3 = S3 / flex_top
    fos_buckling_4 = Tau4 / shear_centroid

    return fos_tension, fos_compression, fos_shear_centroid, fos_shear_glue, fos_buckling_1, fos_buckling_2, fos_buckling_3, fos_buckling_4

def simulate_passes(base_case, step_size):
    lowest_FOS = float('inf')
    train_loads = list(base_case)
    
    # Keep simulating the train passing from left to right until the lowest FOS becomes less than 1
    while lowest_FOS >= 1:
        for i in range(len(train_loads)):
            if i < 2:
                train_loads[i] += 1.35 * 10
            elif i < 4: 
                train_loads[i] += 10
            else:
                train_loads[i] += 1.1 * 10
        print("Simulating load setup: " + str(train_loads))
        print("Total Load: " + str(sum(train_loads)) + " N")
        fos_values = fos(applied_stresses(train_loads, step_size), thin_plate_buckling())
        lowest_FOS = min(fos_values)
        index_FOS = fos_values.index(lowest_FOS)
        
        print("FOS values: " + str(fos_values))
        print("Lowest FOS: " + str(lowest_FOS))
        print("Index of lowest FOS: " + str(index_FOS))
        print("-------------------------")

        
if __name__ == '__main__':
    # Parameter setup
    train_load_loc = (0, 176, 176*1 + 164*1, 176*2 + 164*1, 176*2 + 164*2, 176*3 + 164*2)
    load_setup = (400/6, 400/6, 400/6, 400/6, 400/6, 400/6) # Simplified to have locomotive at the front for every pass
    
    # Global variables
    centroidal_axis = 0
    I = 0
    Q = 0
    Q_glue = 0
    cross_section_height = 0
    centroid_width = 0
    glue_tab_width = 0
    diaphragm_gap = 400
    
    # filename = 'design_zero.txt' 
    # cross_sectional_properties()
    # print("Simulating load setup: " + str(train_loads))
    # print("Total Load: " + str(sum(train_loads)) + " N")
    
    
    # Load setup: Front of train to back of train
    for filename in ['ends.txt', 'sides.txt', 'center.txt']:
        print("\n===== Analyzing " + filename + " =====")
        cross_sectional_properties()
        print("\n----- Starting Simulation -----")
        simulate_passes(load_setup, 1)

    