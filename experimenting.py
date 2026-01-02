"""
Truss in uss
"""
import matplotlib.pyplot as plt

def sfd_envelope(load_setup, step_size):
    # Generating loading positions and loads at each moment
    load_loc_matrix = []
    loads_matrix = []
    front_loc = 0

    while front_loc <= 1200 + 176*3 + 164*2:
        load_loc = [0,0,0,0,0,0]
        loads = [0,0,0,0,0,0]
        for i in range(6):
            if front_loc >= train_load_loc[i] and front_loc - train_load_loc[i] <= 1200:
                load_loc[i] = front_loc - train_load_loc[i]
                loads[i] = load_setup[i]
                
        load_loc_matrix.append(load_loc)
        loads_matrix.append(loads)

        front_loc += step_size
    
    # SFD calculation for each moment
    sfd_matrix = []
    for i in range(len(loads_matrix)):
        rb = 0
        for j in range(6):
            rb += loads_matrix[i][j] * load_loc_matrix[i][j] / 1200
        ra = sum(loads_matrix[i]) - rb
        
        sfd = []
        sfd.append(ra)
        sum_shear = ra
        for j in range(5, -1, -1):
            sum_shear -= loads_matrix[i][j]
            sfd.append(sum_shear)
        
        sfd_matrix.append(sfd)    
    
    plt.figure(figsize=(10, 6))

    # for sfd in sfd_matrix:
    #     plt.plot(x_positions, sfd, linewidth=0.8, alpha=0.35)

    # plt.title("All Shear Force Diagrams (Train Sliding Along Bridge)")
    # plt.xlabel("Position Along Span (mm)")
    # plt.ylabel("Shear Force (kN)")
    # plt.grid(True, linestyle="--", alpha=0.4)

    # plt.show()
    
    
    # # Generating SFD envelope
    # sfd_envelope = []
    # for i in range(1200):
    #     max_shear = -float('inf')
    #     min_shear = float('inf')
    #     for j in range(len(sfd_matrix)):
    #         pos, shear = sfd_matrix[j][i//200], sfd_matrix[j][i//200]
    #         if shear > max_shear:
    #             max_shear = shear
    #         if shear < min_shear:
    #             min_shear = shear
    #     sfd_envelope.append((i, max_shear, min_shear))
    
    # return sfd_envelope
    
    

def bmd_envelope(): 
    pass

if __name__ == '__main__':
    train_load_loc = (0, 176, 176*1 + 164*1, 176*2 + 164*1, 176*2 + 164*2, 176*3 + 164*2)

    sfd_envelope((91, 91, 67.5, 67.5, 67.5, 67.5), 20)