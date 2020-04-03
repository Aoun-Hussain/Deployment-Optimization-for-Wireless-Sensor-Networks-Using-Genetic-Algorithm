import random,math,copy
import networkx as nx;
from networkx.algorithms.approximation.steinertree import steiner_tree

Rad = 15                    ## communication radius of a node
g = 400                     ## cardinality of total sensors in the network
B_SR = 100                  ## Bandwidth between sensor and relay
S_Per_R = 300               ## total sensors per relay
Relay_constraint = 35       ## total relays to be deployed
Sensor_list = []            ## list of sensors
Relay_list = []             ## list of relays
N_R_R_Diction = {}          ## dictionary to show relay-relay connectivity  
N_S_R_Diction = {}          ## dictionary to show sensor-relay connectivity

def main(Population_size, N_GEN, Mutation_prob, Crossover_prob):
    width = 100.0           ## width of the area on which sensors and relays are to be deployed
        
    for i in range(g):      ## generating sensor list
        Sensor_list.append((round(random.uniform(0.0,width),6),round(random.uniform(0.0,width),6)))

    row = 0
    col = 0
    while row <= width:     ## generating relay list
        while col <= width:
            Relay_list.append((float(row),float(col)))
            col += 10
        row += 10
        col = 0

    ###print("Sensor List")
    ###print(Sensor_list)
    ###print("Relay List")
    ###print(Relay_list)
                    
    def distance(a,b):      ## calculates the Euclidean distance
            return math.sqrt((a**2) + (b**2)) 

    def Connection_s_r(SN,RN):  ## calculates the connectivity matrix for sensor x relay layer
            Neighbor_Sensor_Relay = [[[0] for i in range(len(RN))] for j in range(len(SN))]
            for i in range(len(SN)):
                for j in range(len(RN)):
                    dist = distance(abs(SN[i][0] - RN[j][0]), abs(SN[i][1]-RN[j][1]))
                    if dist <= Rad:
                        if i in N_S_R_Diction:
                            N_S_R_Diction[i].append(j)
                        else:
                            N_S_R_Diction[i] = [j]
                        Neighbor_Sensor_Relay[i][j] = 1
                    elif dist > Rad:
                        Neighbor_Sensor_Relay[i][j] = 0
        
            return Neighbor_Sensor_Relay       
            
    N_S_R = (Connection_s_r(Sensor_list,Relay_list))
     
    ###print('Neighbor Matrix of Sensor x Relay Layer', '\n')            
    ###for i in n_s_r:
    ###print(i)
                        
    def Connection_r_r(RN): ## calculates the Connectivity matrix for relay x relay layer
            Neighbor_Relay_Relay = [[[0] for i in range(len(RN))] for j in range(len(RN))]  
            for i in range(len(RN)):
                for j in range(len(RN)):
                    dist = distance(abs(RN[i][0] - RN[j][0]), abs(RN[i][1] - RN[j][1]))
                    if i == j:
                        Neighbor_Relay_Relay[i][j] = 0
                    elif dist <= Rad:
                        if i in N_R_R_Diction:
                            N_R_R_Diction[i].append(j)
                        else:
                            N_R_R_Diction[i] = [j]
                        Neighbor_Relay_Relay[i][j] = 1
                    elif dist > Rad:
                        Neighbor_Relay_Relay[i][j] = 0
            
            return  Neighbor_Relay_Relay       
            
    N_R_R = (Connection_r_r(Relay_list))
     
    ###print('Neighbor Matrix of Relay x Relay Layer', '\n')            
    ###for i in n_r_r:
    ###print(i)

    def Link_s_r(SN,RN): ## calculates the Data link flow matrix for sensor x relay layer
            bandwidth = 100.0
            Linkflow_Sensor_Relay = [[[0] for i in range(len(RN))] for j in range(len(SN))]
            for i in range(len(SN)):
                for j in range(len(RN)):
                    Linkflow_Sensor_Relay[i][j] = bandwidth
            
            return Linkflow_Sensor_Relay

    l_s_r = (Link_s_r(Sensor_list,Relay_list))
        
    ###print('Data Link Flow Matrix of Sensor x Relay Layer', '\n')
    ###for i in l_s_r:
    ###print(i)

    def Link_r_r(RN):   ## calculates the Data link flow matrix for relay x relay layer
            bandwidth = 200.0
            Linkflow_Relay_Relay = [[[0] for i in range(len(RN))] for j in range(len(RN))]
            for i in range(len(RN)):
                for j in range(len(RN)):
                    if i != j:
                        Linkflow_Relay_Relay[i][j] = bandwidth
                    else:
                        Linkflow_Relay_Relay[i][j] = 0
            
            return Linkflow_Relay_Relay

    l_r_r = (Link_r_r(Relay_list))
        
    ###print('Data Link Flow Matrix of Relay x Relay Layer', '\n')
    ###for i in l_r_r:
    ###print(i)
                
    def Energy_s_r(SN,RN):  ## calculates the Energy matrix for sensor x relay layer
            bandwidth = 100.0
            Energy_Sensor_Relay = [[[0] for i in range(len(RN))] for j in range(len(SN))]
            E_radio_S = 50.0 * (10 ** (-9))
            E_radio_R = 100.0 * (10 ** (-9))
            Transmit_amplifier = 100.0 * (10 ** (-12))
            for i in range(len(SN)):
                for j in range(len(RN)):
                    dist = distance(abs(SN[i][0] - RN[j][0]), abs(SN[i][1] - RN[j][1]))
                    energy_sensor_tx = float(bandwidth * (E_radio_S + (Transmit_amplifier * (dist **2))))  ##energy used when sensor transmits data
                    energy_relay_rx = float(bandwidth * E_radio_R)                                         ##energy used when relay receives data 
                    total_energy = energy_sensor_tx + energy_relay_rx
                    Energy_Sensor_Relay[i][j] = total_energy/4                                             ##Only 1/4 sensors will be active per time unit for every relay
                    
            return Energy_Sensor_Relay
        
    e_s_r = (Energy_s_r(Sensor_list,Relay_list))
        
    ###print('Energy Matrix of Sensor x Relay layer', '\n')
    ###for i in e_s_r:
    ###print(i)
    
    def Energy_r_r(RN): ## calculates the Energy matrix for relay x relay layer
            bandwidth = 200.0
            Energy_Relay_Relay = [[[0] for i in range(len(RN))] for j in range(len(RN))]
            E_radio_R = 100.0 * (10 ** (-9))
            Transmit_amplifier = 100.0 * (10 ** (-12))
            for i in range(len(RN)):
                for j in range(len(RN)):
                    dist = distance(abs(RN[i][0] - RN[j][0]), abs(RN[i][1] - RN[j][1]))
                    energy_relay_tx = float(bandwidth * (E_radio_R + (Transmit_amplifier * (dist**2))))  ##energy used when relay transmits data
                    energy_relay_rx = float(bandwidth * E_radio_R)                                       ##energy used when relay receives data 
                    total_energy = (energy_relay_tx + energy_relay_rx)
                    Energy_Relay_Relay[i][j] = total_energy
                    
            return Energy_Relay_Relay

    e_r_r = (Energy_r_r(Relay_list))
        
    ###print('Energy Matrix of Relay x Relay layer', '\n')
    ###for i in e_r_r:
    ###print(i)

    def RelaysDiction(Relay_list):  ## forms the relay graph connectivity its connected edges
        G = nx.Graph()
        for i in range(len(Relay_list)):
            G.add_node(i)
            for j in range(len(Relay_list)):
                if abs(distance(Relay_list[i][0] - Relay_list[j][0],Relay_list[i][1] - Relay_list[j][1])) <= Rad:
                    G.add_edge(i,j)
                    G.add_edge(j,i)
        return G

    RELAYS_GRAPH = RelaysDiction(Relay_list)

    def ProduceOffSprings():
            '''
            this function is used to generate a single individual which contains our system properties i.e. sensor-relay and relay-relay characteristics
            '''
            individual = []                         ## list for an offsprings characteristics
            SR_relation = []                        ## list for the sensr-relay relation of an offsprings characteristics
            RR_relation = []                        ## list for the relay-relay relation of an offsprings characteristics
            relay_record = []                       ## a list to maintain the record of activated relays in the sensr-relay relation of an offsprings characteristics
            
            for i in range(len(Sensor_list)):       ## generating sensor-relay characteristics for an individual 
                lst = [0 for i in range(len(Relay_list))]
                rand = random.choice(N_S_R_Diction[i])
                if (rand) not in relay_record:
                    relay_record.append(rand)
                lst[rand] = 1
                SR_relation.append(lst)
                
            for i in range(len(Relay_list)):        ## generating sensor-relay characteristics for an individual
                lst = [0 for i in range(len(Relay_list))]
                RR_relation.append(lst)                
            for i in relay_record:
                for j in relay_record:
                    if i in N_R_R_Diction[j] and i!=j:
                        RR_relation[i][j] = 1              

            ##            print("SR part of the Individual")
            ##            for j in SR_relation:
            ##                print(j)
            ##            print("RR part of the Individual")
            ##            for j in RR_relation:
            ##                print(j)
            
            individual.append(SR_relation)
            individual.append(RR_relation)
            return individual
        
    def criteria(lst):
            '''
            a function that checks whether the given individual follows the pattern or not.
            '''
            total_linkflow = 0                      ## a variable showing the total linkflow 
            relays_deployed = set()
            lst_SR = lst[0]
            lst_RR = lst[1]
            
            linkflow_list_record = []
            for i in range(len(Relay_list)):
                for j in range(len(Sensor_list)):
                    if lst_SR[j][i] == 1:
                        total_linkflow += l_s_r[j][i]
                        relays_deployed.add(i)
                linkflow_list_record.append(total_linkflow)
                if total_linkflow > B_SR * S_Per_R:
                    return False
                total_linkflow=0
            
                ##max_linkflow = max(linkflow_list_record)
                ##relay_number = linkflow_list_record.index(max_linkflow)
                ##print("Relay Number : ", relay_number, "Max Linkflow : ",max_linkflow,"Sensors connected : ",max_linkflow/100)
            if len(relays_deployed) > Relay_constraint:
                return False
            return True

    def Max_sensors_per_relay(lst):
            '''
            this function tell us the maximum number of sensor per relay.
            '''
            total_linkflow = 0                      ## a variable showing the total linkflow 
            relays_deployed = set()
            lst_SR = lst[0]
            lst_RR = lst[1]
            
            linkflow_list_record = []
            for i in range(len(Relay_list)):
                for j in range(len(Sensor_list)):
                    if lst_SR[j][i] == 1:
                        total_linkflow += l_s_r[j][i]
                linkflow_list_record.append(total_linkflow)
                if total_linkflow > B_SR * S_Per_R:
                    return False
                    break
                total_linkflow=0
            
            max_linkflow = max(linkflow_list_record)
            relay_number = linkflow_list_record.index(max_linkflow)
            print("Maximum number of sensors per relay")
            print("Relay Number : ", relay_number, " Max Linkflow : ",max_linkflow ," Sensors connected : ",max_linkflow/100)
            

    def evaluate(individual):
            '''
            this function checks whether the individual matches the criteria and if yes, then it
            calculates the total energy consumed between sensor-relay and relay-relay.
            '''
            ga_energy = 0
            lst_SR = individual[0]
            lst_RR = individual[1]
            
            if not criteria(individual):            ## condition to check if a individual meets the criteria 
                return float('inf')
            else:
                diction = {}
                for i in range(len(lst_SR)):        ## calculates the energy between sensor-relay communication
                    for j in range(len(lst_SR[0])):
                        if lst_SR[i][j] == 1 and N_S_R[i][j] == 1:
                            ga_energy += e_s_r[i][j]

                for i in range(len(lst_RR)):        ## calculates the energy between relay-relay communication
                    for j in range(len(lst_RR[0])):
                        if lst_RR[i][j] == 1 and N_R_R[i][j] ==1:
                            ga_energy += e_r_r[i][j]/2
                            
                return ga_energy

    def crossover(individual_1,individual_2):
            '''
            the purpose of this function is to crossover between two lists i.e. two individuals
            '''
            rand = random.randint(1,len(individual_1[0])-1)
            relay_record_1 =[]
            relay_record_2 =[]
            
            lst_SR_1 = individual_1[0]
            lst_SR_2 = individual_2[0]
            
            temp_RR_1 = []
            temp_RR_2 = []
            
            temp_SR_1 = copy.deepcopy(lst_SR_1[0:rand] + lst_SR_2[rand:]) 
            temp_SR_2 = copy.deepcopy(lst_SR_2[0:rand] + lst_SR_1[rand:])

            for i in range(len(temp_SR_1)):
                for j in range(len(temp_SR_1[0])):
                    if temp_SR_1[i][j] == 1:
                        relay_record_1.append(j)

            for i in range(len(Relay_list)):
                lst = [0 for i in range(len(Relay_list))]
                temp_RR_1.append(lst)

                
            for i in relay_record_1:
                for j in relay_record_1:
                    if i in N_R_R_Diction[j] and i!=j:
                        temp_RR_1[i][j] = 1

            for i in range(len(temp_SR_2)):
                for j in range(len(temp_SR_2[0])):
                    if temp_SR_2[i][j] == 1:
                        relay_record_2.append(j)

            for i in range(len(Relay_list)):
                lst = [0 for i in range(len(Relay_list))]
                temp_RR_2.append(lst)

                
            for i in relay_record_2:
                for j in relay_record_2:
                    if i in N_R_R_Diction[j] and i!=j:
                        temp_RR_2[i][j] = 1 


            offspring_1 = []
            offspring_2 = []
            
            offspring_1.append(temp_SR_1)
            offspring_1.append(temp_RR_1)
            offspring_2.append(temp_SR_2)
            offspring_2.append(temp_RR_2)
            
            return offspring_1, offspring_2

    def mutate(lst):
            '''
            this function mutates i.e. changes the position of 1 in the given individual.
            '''
            set_of_relays = set()
            Mutated_lst = []
            lst_SR = lst[0]
            lst_RR = lst[1]
            relay_record = []
            breaked = False
            for i in range(len(lst_SR)):                    #Mutating the sensor-relay part of the individual 
                breaked = False
                for j  in N_S_R_Diction[i]:
                    if j in set_of_relays:
                        lst_SR[i][lst_SR[i].index(1)], lst_SR[i][j] = 0,1
                        breaked = True
                        break
                if not breaked:
                    rand = random.choice(N_S_R_Diction[i])
                    lst_SR[i][lst_SR[i].index(1)], lst_SR[i][rand] = 0,1
                    set_of_relays.add(rand)

            for i in range(len(lst_SR)):
                for j in range(len(lst_SR[0])):
                    if lst_SR[i][j] == 1:
                        relay_record.append(j)

            for i in range(len(Relay_list)):
                lst = [0 for i in range(len(Relay_list))]
                lst_RR.append(lst)

                
            for i in relay_record:
                for j in relay_record:
                    if i in N_R_R_Diction[j] and i!=j:
                        lst_RR[i][j] = 1 

            Mutated_lst.append(lst_SR)
            Mutated_lst.append(lst_RR)
            return Mutated_lst

    
    def genetic_algorithm():             
            '''
            this is our primary function which essentially solves our minimum energy consumption problem on UWSNs working on sleep scheduling routing protocols.
            '''
            
            population = [ProduceOffSprings() for i in range(Population_size)]      ## creating individuals for the size of population
            vals = [evaluate(i) for i in population]                                ## evaluating each individual in the population
                
            for i in range(N_GEN):
                vals, population = (zip(*sorted(zip(vals, population))))
                vals = list(vals)
                population = list(population)
            
                for j in range(6):
                    if random.random() < Crossover_prob:
                        population[-j], population[-(j+1)] = crossover(population[j], population[j+1])
                        vals[-j] = evaluate(population[-j])
                        vals[-(j+1)] = evaluate(population[-(j+1)])

                for j in range(6,len(population)-5):
                    if random.random() < Mutation_prob:
                        population[j] = mutate(population[j])
                        vals[j] = evaluate(population[j])
        
            vals, population = (zip(*sorted(zip(vals, population))))
            vals = list(vals)
            population = list(population)
            return population, vals
        
    population, vals = genetic_algorithm()
    relays_deployed = set()

    for i in population:
        ###print("Fittest Individual")
        ###print("SR part of the Individual")
        ###for j in range(len(i[0])):
            ###print(j," : ", i[0][j])
        ###print("RR part of the Individual")
        ###for j in range(len(i[1])):
            ###print(j," : ", i[1][j])
        ###print("Best Individual")
        Max_sensors_per_relay(i)
        sett = set()
        for j in range(len(i[0])):
            sett.add(i[0][j].index(1))
        
        x = steiner_tree(RELAYS_GRAPH, sett)
        #print(len(x.nodes))
        break

    gamma = vals[0]

    Y = 200.0
    E_cosnum_radio_elec_R = 100.0 * 10 ** (-9)
    transmit_amplifier_R = 100.0 * 10 ** (-12)

    for i in range(len(x.nodes) - len(sett)):
        relay_receive_data_rate = Y
        dist = Rad
        energy_relay = (relay_receive_data_rate * E_cosnum_radio_elec_R) + relay_receive_data_rate * (E_cosnum_radio_elec_R+transmit_amplifier_R * (dist ** 2))   
        gamma += energy_relay
    min_energy = gamma 

    print("Maximum relay constraint on the network: ", Relay_constraint)
    print("Maximum candidate location for relays: ", len(Relay_list))
    print("Sensors deployed in the network: ", len(Sensor_list))
    print("Radius of communication of each node: ", Rad)
    print("Relays deployed before steiner nodes: ", len(sett))
    print("Relays deployed after steiner nodes: ", len(x.nodes))
    return (len(x.nodes), min_energy)
                    
nodes = []
energies = []

print("Enter the population size : ")
Population_size = int(input())
print("Enter the number of generations : ")
N_GEN = int(input())
print("Enter the mutation probability : ")
Mutation_prob = float(input())
print("Enter the crossover probability : ")
Crossover_prob = float(input())

for i in range(1):
    node, energy = main(Population_size, N_GEN, Mutation_prob, Crossover_prob)
    nodes.append(node)
    energies.append(energy)


print(sum(nodes)/len(nodes))
print(sum(energies)/len(energies))



    
