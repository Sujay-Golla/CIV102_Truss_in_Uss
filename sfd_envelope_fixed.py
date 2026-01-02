import matplotlib.pyplot as plt

def sfd(load_setup, front_pos):
    # 1. Map train loads to their actual positions on the beam and filter those within [0, 1200]
    loads_on_beam = [] # Stores (position, load_value)
    for i in range(len(load_setup)):
        pos = front_pos - train_load_loc[i]
        if 0 <= pos <= 1200:
            loads_on_beam.append((pos, load_setup[i]))
    
    # Sort loads by their position to ensure correct processing order
    loads_on_beam.sort(key=lambda x: x[0])

    # 2. Calculate reactions Ra and Rb
    sum_moment_about_A = 0 # For calculating Rb (moment about the left support)
    total_loads_sum = 0
    for pos, load in loads_on_beam:
        sum_moment_about_A += load * pos
        total_loads_sum += load
    
    if 1200 == 0: # Avoid division by zero if beam length is 0 (unlikely but good practice)
        rb = 0
    else:
        rb = sum_moment_about_A / 1200
    
    ra = total_loads_sum - rb

    # 3. Construct the Shear Force Diagram for integration (piecewise constant)
    # The format is [(position, shear_value_constant_from_this_position_to_next_change)]
    sfd_for_bmd_points = []
    current_shear = ra # Initial shear at x=0 due to Ra
    sfd_for_bmd_points.append((0, current_shear))
    
    for pos, load in loads_on_beam:
        # Add a point at `pos` with the shear *before* the load (if different from previous point)
        if sfd_for_bmd_points[-1][0] < pos:
            # If there's a segment where shear is constant, add a point at the start of load application.
            sfd_for_bmd_points.append((pos, current_shear))
        
        # Apply the load and update current_shear
        current_shear -= load
        
        # Add a point at `pos` with the shear *after* the load
        # If there was a point at `pos` already (from before load), update it.
        if sfd_for_bmd_points[-1][0] == pos:
            sfd_for_bmd_points[-1] = (pos, current_shear)
        else:
            sfd_for_bmd_points.append((pos, current_shear))
    
    # Ensure the SFD extends to the end of the beam (1200mm) with the final shear value
    # This is critical for correct BMD integration.
    if sfd_for_bmd_points[-1][0] < 1200:
        sfd_for_bmd_points.append((1200, current_shear))
    
    return sfd_for_bmd_points

def sfd_envelope(load_setup, step_size):
    max_shear_envelope = [-float('inf')] * 1201
    min_shear_envelope = [float('inf')] * 1201

    for front_loc in range(0, 1200 + 176*3 + 164*2 + 1, step_size):
        sfd_i = sfd(load_setup, front_loc)
        
        # Interpolate SFD for plotting purposes at `step_size` intervals
        plot_sfd_points = []
        sfd_idx = 0
        current_shear_for_plot = 0
        if sfd_i:
            current_shear_for_plot = sfd_i[0][1] # Shear at x=0

        for x_plot in range(0, 1201):
            while (sfd_idx + 1 < len(sfd_i) and 
                   sfd_i[sfd_idx + 1][0] <= x_plot):
                sfd_idx += 1
                current_shear_for_plot = sfd_i[sfd_idx][1]
            plot_sfd_points.append((x_plot, current_shear_for_plot))

        for pos, shear in plot_sfd_points:
            if 0 <= pos <= 1200:
                idx = int(pos)
                if shear > max_shear_envelope[idx]:
                    max_shear_envelope[idx] = shear
                if shear < min_shear_envelope[idx]:
                    min_shear_envelope[idx] = shear

        plt.plot([pos for pos, shear in plot_sfd_points], [shear for pos, shear in plot_sfd_points], alpha=0.1) # Plot individual SFDs with transparency
        
    plt.title("All Shear Force Diagrams")
    plt.xlabel("Position Along Span (mm)")
    plt.ylabel("Shear Force (N)")
    plt.grid(True, linestyle="--")
    
    plt.plot(range(1201), max_shear_envelope, color='red', linewidth=2, label='Envelope (max)')
    plt.plot(range(1201), min_shear_envelope, color='blue', linewidth=2, label='Envelope (min)')
    plt.xlim(0,1200)
    plt.legend()
    plt.show()

    return max_shear_envelope, min_shear_envelope

def bmd(shear_diagram, step_size):
    bmd_points_critical = [(0, 0)] # Moment at x=0 is 0
    current_moment = 0
    
    # shear_diagram is sorted: [(pos_i, shear_at_pos_i)]
    # shear_at_pos_i is the shear from pos_i to the next pos_i+1
    
    prev_sfd_pos = 0 # Start at 0 for moment calculation
    prev_sfd_shear = 0 # Shear before the first sfd point (effectively 0 for calculating moment from 0)
    
    # Set the initial shear from the first point in `shear_diagram`
    if shear_diagram:
        prev_sfd_shear = shear_diagram[0][1] # This is the shear from 0 to shear_diagram[1][0]
    
    # Integrate each segment of the shear diagram to find moments at critical points
    for i in range(len(shear_diagram)):
        current_sfd_pos = shear_diagram[i][0]
        current_sfd_shear_value = shear_diagram[i][1] # Shear from current_sfd_pos to next point
        
        # Only integrate segments within 0 to 1200
        if prev_sfd_pos < 1200:
            segment_start = prev_sfd_pos
            segment_end = min(current_sfd_pos, 1200) # Integrate up to current_sfd_pos or 1200
            
            if segment_end > segment_start:
                current_moment += prev_sfd_shear * (segment_end - segment_start)
                bmd_points_critical.append((segment_end, current_moment))
        
        prev_sfd_pos = current_sfd_pos
        prev_sfd_shear = current_sfd_shear_value

    # After loop, if prev_sfd_pos < 1200, integrate to 1200 using the last shear.
    if prev_sfd_pos < 1200:
        current_moment += prev_sfd_shear * (1200 - prev_sfd_pos)
        bmd_points_critical.append((1200, current_moment))

    # Now, interpolate to get points at `step_size` intervals for smooth plotting
    final_bmd_plot_points = []
    critical_idx = 0
    
    for x_plot in range(0, 1201, step_size):
        # Find the segment in `bmd_points_critical` that `x_plot` falls into
        while (critical_idx + 1 < len(bmd_points_critical) and
               bmd_points_critical[critical_idx + 1][0] < x_plot):
            critical_idx += 1
        
        if critical_idx < len(bmd_points_critical):
            x1, m1 = bmd_points_critical[critical_idx]
            if critical_idx + 1 < len(bmd_points_critical):
                x2, m2 = bmd_points_critical[critical_idx + 1]
                
                # Linear interpolation
                if x2 > x1:
                    interpolated_moment = m1 + (m2 - m1) * (x_plot - x1) / (x2 - x1)
                else: # x1 == x2 (vertical line in BMD, due to concentrated moment - not in this problem)
                    interpolated_moment = m1 
            else: # Beyond the last critical point, moment is constant at the last value
                interpolated_moment = m1
            
            final_bmd_plot_points.append((x_plot, interpolated_moment))
        else: # No critical points (empty beam, moment is 0)
            final_bmd_plot_points.append((x_plot, 0)) 

    return final_bmd_plot_points


def bmd_envelope(load_setup, step_size):
    max_moment_envelope = [-float('inf')] * 1201
    min_moment_envelope = [float('inf')] * 1201

    # The range for front_loc needs to cover the train moving entirely across the span
    # from its front at 0 to its back leaving the span at 1200.
    # Total length of the train is the last `train_load_loc` value.
    max_front_loc = 1200 + train_load_loc[-1]
    
    for front_loc in range(0, max_front_loc + 1, step_size):
        sfd_for_current_loc = sfd(load_setup, front_loc)
        bmd_for_current_loc = bmd(sfd_for_current_loc, step_size)
            
        # Update envelopes
        for pos_val, moment_val in bmd_for_current_loc:
            # Ensure pos_val is within the span and convert to int for indexing
            if 0 <= pos_val <= 1200:
                idx = int(pos_val)
                if moment_val > max_moment_envelope[idx]:
                    max_moment_envelope[idx] = moment_val
                if moment_val < min_moment_envelope[idx]:
                    min_moment_envelope[idx] = moment_val
            
        # Plot individual BMDs (optional, but requested in original code) with transparency
        plt.plot([pos for pos, moment in bmd_for_current_loc], [moment for pos, moment in bmd_for_current_loc], alpha=0.1)

    plt.title("All Bending Moment Diagrams")
    plt.xlabel("Position Along Span (mm)")
    plt.ylabel("Moment (Nmm)")
    plt.grid(True, linestyle="--")
    
    # Plot the envelopes
    plt.plot(range(1201), max_moment_envelope, color='red', linewidth=2, label='Envelope (max)')
    plt.plot(range(1201), min_moment_envelope, color='blue', linewidth=2, label='Envelope (min)')
    plt.xlim(0,1200)
    plt.legend()
    plt.show()

    return max_moment_envelope, min_moment_envelope

if __name__ == '__main__':
    train_load_loc = (0, 176, 176*1 + 164*1, 176*2 + 164*1, 176*2 + 164*2, 176*3 + 164*2)
    # Load setup: Front of train to back of train

    sfd_envelope((91, 91, 67.5, 67.5, 67.5, 67.5), 1)
    bmd_envelope((91, 91, 67.5, 67.5, 67.5, 67.5), 1)