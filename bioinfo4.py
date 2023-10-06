import numpy as np
import copy
#WEEK1
def is_an_edge(node_1, node_2):
    edge = node_1 + "->" + node_2
    if edge in list_of_edges:
        return True
    else:
        return False

def is_leave(node):
    return node in list_of_leaves


def find_distance(node, visited=[],length=0):
    # Bad function
    count =0 
    count +=1
    visited.append(node)
    for node_2 in list_of_edges[node]:
        if node_2 in visited:
            continue
        else:
            edge = node + "->" + node_2            
            path = leave + "->" + node_2
            path_reverse = node_2 + "->" + leave
            if is_leave(node_2) and path not in list_of_paths:
                list_of_paths[path] = length  + list_of_lengths[edge]
                list_of_paths[path_reverse] = length + list_of_lengths[edge]
                continue
        find_distance(node_2,visited,length + list_of_lengths[edge])
    if count > 9999:
        return None
    return None


#Find the parental node of leaves and its length
def LimbLengthList(D, j):
    # Create a list of length
    assert D.shape[0] == D.shape[1], "Matrix doesn't have the same numbers of rows and columns"
    n = D.shape[0]
    list_of_lengths = dict()
    for i in range(n):
        if i == j:
            continue
        for k in range(n):
            if k == j or k == i:
                continue
            limb = (D[i,j] + D[j,k] - D[i,k])/2
            list_of_lengths[f"{i}-{k}"]=limb
    return list_of_lengths
def FindNodesNearestToJ(D,j):
    # Find the nearest node to leave
    list_of_lengths = LimbLengthList(D,j)
    return min(list_of_lengths.keys(), key=(lambda k: abs(list_of_lengths[k])))
def LimbLength(D,j):
    #Find the length of the nearest node to leave
    list_of_lengths = LimbLengthList(D,j)
    return min(list_of_lengths.values(), key=(lambda k: abs(k)))





def TString2Dict(T_str):
    #Transform input of tree from string to dictionary
    T = dict()
    T_str = T_str.strip().split("\n")
    for i in T_str:
        edge,length = i.split(":")
        T[edge] = int(length)
    return T

def T2EdgesList(T):
    #Transform the tree to list of edges
    if type(T) is str:
        T = TString2Dict(T)
    list_of_edges = {}
    for edge in T:
        node1,node2 = edge.split("->")
        if node1 in list_of_edges:
            list_of_edges[node1].append(node2)
        else:
            list_of_edges[node1] = [node2]
    return list_of_edges

def PathOf2Leaves(list_of_edges,leaf_start,leaf_end,visited = None):
    #Find the unique path between 2 leaves
    if visited == None:
        visited = []
    visited.append(leaf_start)
    if leaf_end in list_of_edges[leaf_start]:
        return [leaf_end,leaf_start]
    
    for node in list_of_edges[leaf_start]:
        if node in visited:
            continue
        P = PathOf2Leaves(list_of_edges,node,leaf_end,visited)
        if P == None:
            continue
        else:
            P.append(leaf_start)
        return P



def NodeExits(T,i,k, D_ix):
    #Find if the node X with distancec D_ix from i alrady exist in the path between i and k
    #if node exist, return its name
    #if not, return the position of the new node under format (node1,node2,distance_node1_to_X)
    i = str(i)
    k = str(k)
    list_of_edges = T2EdgesList(T)
    path_ik = PathOf2Leaves(list_of_edges,i,k)
    path_ik = path_ik[::-1]
    D_node_to_i = 0
    for cursor, node in enumerate(path_ik[:-1]):
        current_edge = node+"->"+path_ik[cursor+1]
        D_node_to_i += T[current_edge]
        if node == i:
            if D_node_to_i == D_ix:
                return True,path_ik[cursor+1]
            elif D_node_to_i > D_ix:
                return False,path_ik[cursor],path_ik[cursor+1],D_ix
            else:
                continue
        else:
            if D_node_to_i == D_ix:
                return True,path_ik[cursor+1]
            elif D_node_to_i > D_ix:
                return False,path_ik[cursor], path_ik[cursor+1], T[current_edge] - (D_node_to_i - D_ix)
            else:
                continue 
            
    return False,None


def MkEdge(node1,node2):
    return f"{node1}->{node2}"

def Make2NodesTree(D):
    assert D.shape  == (2,2), f"Matrix does not have the required shape of (2,2) but {D.shape} "
    T = dict()
    n = D.shape[0]
    for i in range(n):
        for k in range(n):
            if i != k:
                T[MkEdge(i,k)] = D[i,k]
    return T

def RvEdge(edge):
    node1,node2 = edge.split("->")
    return node2 + "->" + node1

def AddEdge2T(T,edge,length):
    rv_edge = RvEdge(edge)
    T[edge] = length
    T[rv_edge] = length
    return T

def AddNode2T(T,i,k,d_i2node,name=None):
    T = T.copy()
    edge_to_split = MkEdge(i,k)
    reverse_edge_to_split = RvEdge(edge_to_split)
    d_node2k = T[edge_to_split] - d_i2node
    if name == None:
        new_node = "n" + i + k
    else:
        new_node = name
    new_edge1 = MkEdge(i,new_node)
    new_edge2 = MkEdge(new_node,k)
    T = AddEdge2T(T,new_edge1,d_i2node)
    T = AddEdge2T(T,new_edge2,d_node2k)
    del T[edge_to_split]
    del T[reverse_edge_to_split]
    return T

def SortTree(T):
    T = {k:v for k,v in sorted(T.items(),key=(lambda item: item[0]))}
    T = {k:v for k,v in sorted(T.items(),key=(lambda item: int(item[0].split("->")[0])))}
    return T


def AdditivePhylogeny(D, rc_count= None):
    #Recursive algorithm to return a phylogenic tree from a Distance matrix

    ###recursive count and new node count 
    if rc_count == None:
        rc_count = 0
    else:
        rc_count +=1
    global nn_count
    if rc_count == 0:
        nn_count = 0

    ### Assure that distance matrix is symmetrical
    assert D.shape[0] == D.shape[1], "Matrix doesn't have the same numbers of rows and columns"

    n = D.shape[0]
    ### The base
    if n == 2:
        return Make2NodesTree(D)
    
    limbLength = LimbLength(D,n-1) #0-based range
    for j in range(n-1): # To avoid negative result since if j = n then D[n,j] == D[j,n] == 0 
        D[j,n-1] = D[j,n-1] - limbLength
        D[n-1,j] = D[j,n-1]
    (i,k) = FindNodesNearestToJ(D,n-1).split("-")
    x = D[int(i),n-1]
    D = np.delete(np.delete(D,n-1,0),n-1,1)
    T = AdditivePhylogeny(D,rc_count)
    int_node_exist = NodeExits(T,i,k,x)
    
    if int_node_exist[0] == True:
        int_node = int_node_exist[1]
        new_limb = int_node+"->"+str(n-1)
        reverse_new_limb = RvEdge(new_limb)
        T[new_limb] = limbLength
        T[reverse_new_limb] = limbLength
    else:
        node1,node2,D_nn2n1 = int_node_exist[1:]
        name = str(n+rc_count+nn_count)
        nn_count += 1
        T = AddNode2T(T,node1,node2,D_nn2n1,name)
#         new_int_node = "n" + node1+node2
        new_int_node = name
        new_limb = new_int_node+"->"+str(n-1)
        reverse_new_limb = RvEdge(new_limb)
        T[new_limb] = limbLength
        T[reverse_new_limb] = limbLength
    return T


####WEEK 2 #########
def input_to_D_array(ip):
    ip = ip.strip()
    ip = ip.split('\n')
    n = int(ip.pop(0))
    D_array = np.zeros((n,n))
    try:
        for i,row in enumerate(ip):
            D_array[i] = [int(k) for k in row.strip().split(" ")]
    except:
        for i,row in enumerate(ip):
            D_array[i] = [int(k) for k in row.strip().split("\t")]
    return D_array

def D_array_to_dict(D_array):
    n = D_array.shape[0]
    D_dict = {}
    for i in range(n-1):
        for k in range(i+1, n):
            D_dict[f"{i}-{k}"] = D_array[i,k]
    return D_dict

def input_to_D_dict(ip):
    D_array = input_to_D_array(ip)
    return D_array_to_dict(D_array)

def find_closest_clusters(D_dict):
    return min(D_dict.keys(), key=(lambda k: D_dict[k]))

def UPGMA(D_array,n,nn_count =None):
    D_dict = D_array_to_dict(D_array)

    Clusters = {}
    for i in range(n):
        Clusters[str(i)] = [str(i)]
    Clusters_ref = Clusters.copy()

    T = {}

    Age = {}
    for node in Clusters:
        Age[node] = 0
    if nn_count == None:
        nn_count = 1
    while len(Clusters) > 2 :

        # find closest clusters
        closest_clusters = find_closest_clusters(D_dict)
        cluster_1, cluster_2 = closest_clusters.split("-")

        if closest_clusters in D_dict:
            closest_distance = D_dict[closest_clusters]
        #update Distance and Clusters
        del D_dict[f"{cluster_1}-{cluster_2}"], Clusters[cluster_1], Clusters[cluster_2]

        #new_node = f"{cluster_1},{cluster_2}"
        new_node = str(n - 1 + nn_count)
        nn_count += 1
        Clusters_ref[new_node] = Clusters_ref[cluster_1] + Clusters_ref[cluster_2]    
        #Merge 2 cluster and update distance
        for i in Clusters:
            if f"{i}-{cluster_1}" in D_dict:
                Di_to_cluster1 = D_dict[f"{i}-{cluster_1}"].copy()
                del D_dict[f"{i}-{cluster_1}"]
            else:
                Di_to_cluster1 = D_dict[f"{cluster_1}-{i}"]
                del D_dict[f"{cluster_1}-{i}"]

            if f"{i}-{cluster_2}" in D_dict:
                Di_to_cluster2 = D_dict[f"{i}-{cluster_2}"].copy()
                del D_dict[f"{i}-{cluster_2}"]
            else:
                Di_to_cluster2 = D_dict[f"{cluster_2}-{i}"]
                del D_dict[f"{cluster_2}-{i}"]

            Di_new_node = (
                            (Di_to_cluster1*len(Clusters_ref[cluster_1]))+
                           (Di_to_cluster2*len(Clusters_ref[cluster_2]))
                           )/(len(Clusters_ref[cluster_1])+len(Clusters_ref[cluster_2]))    

            D_dict[f"{i}-{new_node}"] = Di_new_node

        Clusters[new_node] = Clusters_ref[cluster_1] + Clusters_ref[cluster_2]

        # update age, add new node to T and connect new node 
        Age[new_node] = closest_distance/2
        new_edge = f"{new_node}->{cluster_1}"
        edge_length = Age[new_node]- Age[cluster_1]
        T = AddEdge2T(T,new_edge,edge_length)

        new_edge = f"{new_node}->{cluster_2}"
        edge_length = Age[new_node]- Age[cluster_2]
        T = AddEdge2T(T,new_edge,edge_length)



    root =str(n - 1 + nn_count)
    closest_clusters = find_closest_clusters(D_dict)
    closest_distance = D_dict[closest_clusters]
    cluster_1, cluster_2 = closest_clusters.split("-")

    Age[root] = D_dict[find_closest_clusters(D_dict)]/2
    new_edge = f"{root}->{cluster_1}"
    edge_length = Age[root]- Age[cluster_1]
    T = AddEdge2T(T,new_edge,edge_length)

    new_edge = f"{root}->{cluster_2}"
    edge_length = Age[root]- Age[cluster_2]
    T = AddEdge2T(T,new_edge,edge_length)
    T = SortTree(T)
    return T

def FindTotalDistance(D_dict):
    Total_distance = {}
    for key in D_dict:
        node1, node2 = key.split("-")
        if node1 in Total_distance :
            Total_distance[node1] = Total_distance[node1] + D_dict[key]
        else:
            Total_distance[node1] = D_dict[key]

        if node2 in Total_distance :
            Total_distance[node2] = Total_distance[node2] + D_dict[key]
        else:
            Total_distance[node2] = D_dict[key]
    return Total_distance

def FindN(D_dict):
    x = len(D_dict)
    count = 1
    while x > 0:
        x -= count
        count += 1
    return count

def NJD(D_dict, Total_distance = None):
    if Total_distance == None:
        Total_distance = FindTotalDistance(D_dict)
    D_copy = D_dict.copy()
    n = FindN(D_dict)
    for key in D_copy:
        node1,node2 = key.split("-")
        new_distance = (n-2)*D_dict[key] - Total_distance[node1] - Total_distance[node2]
        D_copy[key] = new_distance
    return D_copy

def RemoveNode(D_dict,node):
    D_copy = copy.deepcopy(D_dict)
    keys = list(D_copy.keys())
    for key in keys:
        if node in key.split("-"):
            del D_copy[key]
    return D_copy

def FindDistance(D_dict,i,j):
    try:
        return D_dict[f"{i}-{j}"]
    except:
        return D_dict[f"{j}-{i}"]
    


def NeighborJoining(D_dict, rc_count=None):
    D_dict = copy.deepcopy(D_dict)
    ###recursive count and new node count 

    if rc_count == None:
        rc_count = 0
        global nn_count
        nn_count = 0

    else:
        rc_count +=1
    n = FindN(D_dict)
    #base
    if len(D_dict) == 1:
        T = {}
        #Make 2-node Tree
        node1, node2 = list(D_dict.keys())[0].split("-")
        edge = MkEdge(node1,node2)
        length = list(D_dict.values())[0]
        length = round(length,3)
        rv_edge = RvEdge(edge)
        T[edge] = length
        T[rv_edge] = length
        return T
    Total_distance = FindTotalDistance(D_dict)
    NJD_dict = NJD(D_dict, Total_distance)
    neighbors_ij = min(NJD_dict.keys(),key=(lambda k: NJD_dict[k]))
    node_i, node_j = neighbors_ij.split("-")
    delta_ij = (Total_distance[node_i] - Total_distance[node_j])/(n-2)
    limbLength_i = round((D_dict[neighbors_ij] + delta_ij)/2,3)
    limbLength_j = round((D_dict[neighbors_ij] - delta_ij)/2,3)
    #
    new_node = str(n + rc_count + nn_count)
    nn_count += 1
    for node_k in Total_distance:
        if node_k != node_i and node_k != node_j:
            D_i_k = FindDistance(D_dict, node_i, node_k)
            D_j_k = FindDistance(D_dict, node_j, node_k)
            D_nn = D_dict[neighbors_ij]
            D_k_nn = (D_i_k + D_j_k - D_nn)/2
            D_dict[f"{node_k}-{new_node}"] = D_k_nn
    #
    D_dict = RemoveNode(D_dict,node_i)
    D_dict = RemoveNode(D_dict,node_j)
    
    T = NeighborJoining(D_dict,rc_count)
    
    #
    edge1 = MkEdge(new_node, node_i)
    rv_edge1 = RvEdge(edge1)
    T[edge1] = limbLength_i
    T[rv_edge1] = limbLength_i

    edge2 = MkEdge(new_node, node_j)
    rv_edge2 = RvEdge(edge2)
    T[edge2] = limbLength_j
    T[rv_edge2] = limbLength_j
    
    return T


####WEEEK 3
def HammingDistance(p,q):
    return sum([(i!=k) for i,k in zip(p,q)])
def IsAlphabet(string):
    try:
        int(string)
    except:
        return True
    
def ParsingInputW3(ip, i = None):
    T = {}
    T["nodes"] = {}
    T["characters"] = {}
    T["adjacency_list"] = {}
    T["leaves"] = {}
    ip = ip.strip()
    ip = ip.split("\n")
    n = ip.pop(0)
    for edge in ip:
        node1,node2 = edge.split("->")
        if i != None:
            if IsAlphabet(node2):
                node2 = node2[:i]
        
        if IsAlphabet(node2):
            leave_counts = len(T["leaves"])
            leave_counts = str(leave_counts)
            T["leaves"][leave_counts] = [node2]
            T["characters"][leave_counts] = node2
            T["nodes"][leave_counts] = []
            if node1 not in T["nodes"]:
                T["nodes"][node1] = [leave_counts]
                T["characters"][node1] = ""
            else:
                T["nodes"][node1].append(leave_counts)
        else:
            if node1 not in T["nodes"]:
                T["nodes"][node1] = [node2]
                T["characters"][node1] = ""

            else:
                T["nodes"][node1].append(node2)
    T["nodes"] = {k:v for k,v in sorted(T["nodes"].items(), key=(lambda item: int(item[0])))}
    T["characters"] = {k:v for k,v in sorted(T["characters"].items(), key=(lambda item: int(item[0])))}    
    return T

def IsRipe(node_v,T,Tag):
    s = map(lambda child_node: Tag[child_node], T["nodes"][node_v])
    if Tag[node_v] == 0 and sum(s) == 2 :
        return True
    else:
        return False

    
def UntagExists(Tag):
    if len(Tag) > sum(Tag.values()):
        return True
    else:
        return False

def FindRipeList(T,Tag):
    ripe_list = []
    if UntagExists(Tag):
        for node_v in Tag:
            if Tag[node_v] == 0 and IsRipe(node_v,T,Tag):
                ripe_list.append(node_v)
    return ripe_list


def SmallParsimony(T, Character):
    alphabet = ["A","C","G","T"]
    Tag = {}
    S = {}
    for node_v in T["nodes"]:
        Tag[node_v] = 0
        S[node_v] = {}
        if node_v in T["leaves"]:
            Tag[node_v] = 1
            for k in alphabet:
                if Character[node_v] == k:
                    S[node_v][k] = 0
                else:
                    S[node_v][k] = float("inf")
    Backtrack = {}
    RipeList = FindRipeList(T,Tag)
    while len(RipeList) > 0:
        node_v = RipeList[0]
        son,daughter = T["nodes"][node_v]
        Tag[node_v] = 1
        for k in alphabet:
            son_min_value = np.min([S[son][j] + HammingDistance(j,k) for j in alphabet])
            son_min_char_idx = np.argmin([S[son][j] + HammingDistance(j,k) for j in alphabet])
            son_min_char = alphabet[son_min_char_idx]

            daughter_min_value = np.min([S[daughter][i] + HammingDistance(i,k) for i in alphabet])
            daughter_min_char_idx = np.argmin([S[daughter][i] + HammingDistance(i,k) for i in alphabet])
            daughter_min_char = alphabet[daughter_min_char_idx]

            S[node_v][k] = son_min_value+daughter_min_value
            


            
        #Update RipeList
        RipeList = FindRipeList(T,Tag)

    return S

def BackTrack(T,S):
    alphabet = ["A","C","G","T"]
    Backtrack = {}
    score = 0
    for node in range(len(T["nodes"])-1,len(T["leaves"])-1,-1):
        if node == len(T["nodes"])-1:
            root = str(node)
            root_node_char = min(S[root],key=S[root].get)
            score = S[root][root_node_char]
            Backtrack[root] = root_node_char

            son,daughter = T["nodes"][root]        
            k = root_node_char
            son_min_value = np.min([S[son][j] + HammingDistance(j,k) for j in alphabet])
            son_min_char_idx = np.argmin([S[son][j] + HammingDistance(j,k)*1.01 for j in alphabet])
            son_min_char = alphabet[son_min_char_idx]

            daughter_min_value = np.min([S[daughter][i] + HammingDistance(i,k) for i in alphabet])
            daughter_min_char_idx = np.argmin([S[daughter][i] + HammingDistance(i,k)*1.01 for i in alphabet])
            daughter_min_char = alphabet[daughter_min_char_idx]
            Backtrack[son] = son_min_char
            Backtrack[daughter] = daughter_min_char
        else:
            node = str(node)
            son,daughter = T["nodes"][node]        
            k = Backtrack[node]
            son_min_value = np.min([S[son][j] + HammingDistance(j,k) for j in alphabet])
            son_min_char_idx = np.argmin([S[son][j] + HammingDistance(j,k)*1.01 for j in alphabet])
            son_min_char = alphabet[son_min_char_idx]

            daughter_min_value = np.min([S[daughter][i] + HammingDistance(i,k) for i in alphabet])
            daughter_min_char_idx = np.argmin([S[daughter][i] + HammingDistance(i,k)*1.01 for i in alphabet])
            daughter_min_char = alphabet[daughter_min_char_idx]
            Backtrack[son] = son_min_char
            Backtrack[daughter] = daughter_min_char
    return score,Backtrack

def SmallParsimonyString(T):
    
    P_score = 0

    for i_character in range(len(T["leaves"]["0"][0])):
        Character = {}
        for (key, value) in T["leaves"].items():
            Character[key] = value[0][i_character]

        S = SmallParsimony(T,Character)
        score,Backtrack = BackTrack(T,S)
        for node in T["nodes"]:
            if node not in T["leaves"]:
                T["characters"][node] += Backtrack[node]
        P_score += score
    T["P_score"] = P_score
    for node_v in T["nodes"]:
        if node_v not in T["leaves"]:
            for node_child in T["nodes"][node_v]:
                node_v_char = T["characters"][node_v]
                node_child_char = T["characters"][node_child]
                edge = MkEdge(node_v_char,node_child_char)
                rv_edge = RvEdge(edge)
                length = HammingDistance(node_v_char, node_child_char)
                T["adjacency_list"][edge] = length
                T["adjacency_list"][rv_edge] = length
                
    return T

def TestTree(T):
    sum_HD = 0
    for k,v in T["adjacency_list"].items():
        sum_HD +=v
    if int(sum_HD/2) == T["P_score"]:
        print ("PASSED TEST")
    else:
        print ("DID NOT PASS")


def ParsingInputW3_Unrooted(ip, i = None):
    T = {}
    T["nodes"] = {}
    T["characters"] = {}
    T["adjacency_list"] = {}
    T["leaves"] = {}
    T["edges"] = {}
    nodes = []
    ip = ip.strip().split("\n")
    n = int(ip.pop(0))

    for edge in ip:
        node1,node2 = edge.split("->")

        if node1 not in nodes:
            nodes.append(node1)
        if node2 not in nodes:
            nodes.append(node2)

    for node in nodes :
        if IsAlphabet(node):
            leave_counts = len(T["leaves"])
            leave_counts = str(leave_counts)
            T["leaves"][node] = leave_counts
            T["characters"][leave_counts] = node
            T["nodes"][leave_counts] = []
        else:
            T["characters"][node] = ""

    for edge in ip:
        node1,node2 = [T["leaves"][node] if IsAlphabet(node) else node for node in edge.split("->")]
        if node1 not in T["edges"]:
            T["edges"][node1] = [node2]
        else:
            T["edges"][node1].append(node2)
    T["leaves"] = {v:[k] for k,v in T["leaves"].items()}

    return T

def AddRoot(T):
    #Add root
    new_root = str(len(T["edges"]))
    son = str(len(T["edges"])-1)
    daughter = str(len(T["edges"])-2)
    assert son in T["edges"][daughter]
    #Add new edges, Delete old edges
    T["characters"][new_root] = ""
    T["edges"][new_root]= [son,daughter]
    T["edges"][son].remove(daughter)
    T["edges"][daughter].remove(son)
    
    #Log of change
    T["changes"] = {"new_root":new_root,"son":son,"daughter":daughter}
    
    #Delete all backward edges
    z = [son,daughter]
    count = 0
    while len(z) > 0 :
        count += 1
        new_z = []
        for i in z:
            for ii in T["edges"][i]:
                T["edges"][ii].remove(i)
            new_z +=  T["edges"][i]
        z = new_z
        if count > 9999999:
            print("Time Out", count)
            break
    #polishing        
    T["nodes"] = {k:v for k,v in sorted(T["edges"].items(), key=(lambda item: int(item[0])))}
    T["characters"] = {k:v for k,v in sorted(T["characters"].items(), key=(lambda item: int(item[0])))}
    return T

def ReverseRootedTree(T):
    new_root = T["changes"]["new_root"]
    son = T["changes"]["son"]
    daughter =T["changes"]["daughter"]

    del T["nodes"][new_root]
    T["nodes"][son].append(daughter)
    del T["changes"]
    return T

def SmallParsimonyString_Unrooted(T):
    T = AddRoot(T)
    P_score = 0

    for i_character in range(len(T["leaves"]["0"][0])):
        Character = {}
        for (key, value) in T["leaves"].items():
            Character[key] = value[0][i_character]

        S = SmallParsimony(T,Character)
        score,Backtrack = BackTrack(T,S)
        for node in T["nodes"]:
            if node not in T["leaves"]:
                T["characters"][node] += Backtrack[node]
        P_score += score
    T["P_score"] = P_score
    T = ReverseRootedTree(T)
    for node_v in T["nodes"]:
        if node_v not in T["leaves"]:
            for node_child in T["nodes"][node_v]:
                node_v_char = T["characters"][node_v]
                node_child_char = T["characters"][node_child]
                edge = MkEdge(node_v_char,node_child_char)
                rv_edge = RvEdge(edge)
                length = HammingDistance(node_v_char, node_child_char)
                T["adjacency_list"][edge] = length
                T["adjacency_list"][rv_edge] = length
    return T


#######WEEk4#####
aa_list = """G 57
A 71
S 87
P 97
V 99
T 101
C 103
I 113
L 113
N 114
D 115
K 128
Q 128
E 129
M 131
H 137
F 147
R 156
Y 163
W 186"""
aa_list = aa_list.strip().split("\n")
aa_weights = {}
for line in aa_list:
    aa,w = line.split(" ")
    aa_weights[aa] = int(w)
aa_weights_rv = {v:k for k,v in aa_weights.items()}

def Graph(full_spectrum):
    #Full spectrum contain zero and the mass of the whole peptide
    if 0 not in full_spectrum:
        full_spectrum.append(0)
    full_spectrum = sorted(full_spectrum)
        
    n = len(full_spectrum)
    arr = np.zeros([n,n])
    for i in range(len(full_spectrum)-1):
        for k in range(i,len(full_spectrum)):
            arr[i,k]=abs(full_spectrum[i]-full_spectrum[k])
    graph = {}
    for w in aa_weights_rv:
        result = np.where(arr==w)
        result = list(zip(result[0],result[1]))
        for x,y in result:
            graph[f"{full_spectrum[x]}->{full_spectrum[y]}"] = aa_weights_rv[w]
    graph =  {k:v for k,v in sorted(graph.items(),key=(lambda item: int(item[0].split("->")[0])))}
    return graph

def Adjency(full_spectrum):
    #Full spectrum contain zero and the mass of the whole peptide    
    if 0 not in full_spectrum:
        full_spectrum.append(0)
    full_spectrum = sorted(full_spectrum)
            
    n = len(full_spectrum)
    arr = np.zeros([n,n])
    for i in range(len(full_spectrum)-1):
        for k in range(i,len(full_spectrum)):
            arr[i,k]=abs(full_spectrum[i]-full_spectrum[k])
    adj = {}
    for w in aa_weights_rv:
        result = np.where(arr==w)
        result = list(zip(result[0],result[1]))
        for x,y in result:
            if x not in adj:
                adj[x] = [y]
            else:
                adj[x].append(y)
    adj =  {k:v for k,v in sorted(adj.items(),key=(lambda item: int(item[0])))}
    return adj


def FindAllPaths(start_node,end_node,adjency):
    adjency = copy.deepcopy(adjency)
    succeeded = {}
    visited = []
    all_paths = []
    path = []

    current_node = start_node


    while True :
        path.append(current_node)
        #if reach destination, save path and dynamic programming to save possible paths for each node
        if current_node == end_node:            
            all_paths.append(path)
            for i in range(len(path) -1):
                node_i = path[i]
                path_i = path[i+1:]
                if node_i not in succeeded:
                    succeeded[node_i] = [path_i]
                elif path_i not in succeeded[node_i]:
                    succeeded[node_i].append(path_i)
            current_node = path[-2]
            path = path[:-2]
            continue
        #if reach already visited node and this node doesn't lead to destination, move back
        if current_node in visited and current_node not in succeeded:
            visited.append(current_node)
            current_node = path[-2]
            path = path[:-2]
            continue

            
        #if reach a visited node and the path from this node has already been saved, populate new path from this save path 
        if current_node in visited and current_node in succeeded:
            for succeeded_path in succeeded[current_node]:
                new_path = path + succeeded_path
                all_paths.append(new_path)
                for i in range(len(path) -1):
                    node_i = path[i]
                    path_i = path[i+1:] + succeeded_path
                    if node_i not in succeeded:
                        succeeded[node_i] = [path_i]
                    elif path_i not in succeeded[node_i]:
                        succeeded[node_i].append(path_i)
            current_node = path[-2]
            path = path[:-2]
            continue
            
        #if reach dead end, move back a node, if dead end is start node, return all_paths
        if current_node not in adjency or len(adjency[current_node]) == 0:
            visited.append(current_node)
            if current_node == start_node:
                return all_paths
            current_node = path[-2]
            path = path[:-2]
            continue
            
        #move forward, pop a node in adjency list
        if current_node in adjency and len(adjency[current_node]) > 0:
            current_node = adjency[current_node].pop(-1)
            continue

def Weight(peptide):
    w = 0
    for aa in peptide:
        w += aa_weights[aa]
    return w 

def IdealSpectrum(peptide):
    prefix_peptides = []
    suffix_peptides = []
    for i in range(1,len(peptide)):
        prefix_peptides.append(peptide[:i])
        suffix_peptides.append(peptide[i:])
    fragments = prefix_peptides + suffix_peptides
    spectrum = [Weight(fragment) for fragment in fragments]
    spectrum.append(Weight(peptide))
    spectrum.append(0)
    return sorted(spectrum)


def DecodeIdealSpectrum(full_spectrum):
    print("Make sure the full spectrum is entered")
    if 0 not in full_spectrum:
        full_spectrum.append(0)
    full_spectrum = sorted(full_spectrum)
    adj = Adjency(full_spectrum)
    nodes = []
    for k,v in adj.items():
        nodes.append(k)
        nodes += v
    nodes = list(set(nodes))
    start = nodes[0]
    end = nodes[-1]
    print(start,end)
    all_paths = FindAllPaths(start,end,adj)

    graph = Graph(full_spectrum)

    ideal_peptides = []
    for path in all_paths:
        peptide = ""
        for i in range(len(path)-1):
            prefix = full_spectrum[path[i]]
            suffix = full_spectrum[path[i+1]]
            edge = f"{prefix}->{suffix}"
            peptide += graph[edge]
        if IdealSpectrum(peptide) == full_spectrum:
            ideal_peptides.append(peptide)
    return ideal_peptides


def PepToVector(peptide):
    prefix_peptides = []
    for i in range(1,len(peptide)+1):
        prefix_peptides.append(peptide[:i])
    prefix_masses = [Weight(p) for p in prefix_peptides]
    binary_mass_vector = np.zeros(prefix_masses[-1])
    for i in prefix_masses:
        binary_mass_vector[i-1] = 1
    return binary_mass_vector

def VectorToPep(binary_mass_vector):
    print("Make sure to add 0 in the first position in input")
    result = np.where(binary_mass_vector==1)
    spectrum = [ mass for mass in result[0]]
    suffix_mass = []
    for mass in spectrum[:-1]:
        suffix_mass.append(spectrum[-1] - mass)
    spectrum += suffix_mass
    peptides = DecodeIdealSpectrum(spectrum)
    return peptides


def FindAllPathWeightedDAG(start_node,end_node,adjency,real_spectra,poison = None):
    if poison == None:
        poison = []
        
    succeeded = {}
    visited = []
    all_paths = []
    path = []

    current_node = start_node
    while True :
        path.append(current_node)
        #if reach destination, save path and calculate weight(score) of the path for each node in the path
        #if node already has a path, chose the path with higher weight
        if current_node == end_node:
            all_paths.append(path)
            for i in range(len(path) -1):
                node_i = path[i]
                path_i = path[i+1:]
                if node_i not in succeeded:
                    new_sum = sum([real_spectra[k] for k in path_i])
                    succeeded[node_i] = [path_i,new_sum]
                else:
                    new_sum = sum([real_spectra[k] for k in path_i])
                    old_sum = succeeded[node_i][1]
                    if new_sum > old_sum:
                        succeeded[node_i] = [path_i,new_sum]
            current_node = path[-2]
            path = path[:-2]
            continue
            
        #if reach already visited node and this node doesn't lead to destination, move back
        if current_node in visited and current_node not in succeeded:
            visited.append(current_node)
            current_node = path[-2]
            path = path[:-2]
            continue
            
        #if reach a poison node (retristed), move back, poison node become visited
        if current_node in poison:
            visited.append(current_node)
            current_node = path[-2]
            path = path[:-2]
            continue
            
        #if reach a visited node and the path from this node has already been saved, populate new path from this save path
        #calculate weight(score) of the path for potentially new node in the path
        #if node already has a path, chose the path with higher weight        
        if current_node in visited and current_node in succeeded:
            succeeded_path = succeeded[current_node][0]
            new_path = path + succeeded_path
            all_paths.append(new_path)
            for i in range(len(path) -1):
                node_i = path[i]
                path_i = path[i+1:] + succeeded_path
                if node_i not in succeeded:
                    new_sum = sum([real_spectra[k] for k in path_i])
                    succeeded[node_i] = [path_i,new_sum]
                else:
                    new_sum = sum([real_spectra[k] for k in path_i])
                    old_sum = succeeded[node_i][1]
                    if new_sum > old_sum:
                        succeeded[node_i] = [path_i,new_sum]

            current_node = path[-2]
            path = path[:-2]
            continue
       
        #if reach dead end, move back a node, if dead end is start node, return all_paths
        if current_node not in adjency or len(adjency[current_node]) == 0:
            visited.append(current_node)
            if current_node == start_node:
                return all_paths
            current_node = path[-2]
            path = path[:-2]
            continue
        
        #move forward, pop a node in adjency list
        if current_node in adjency and len(adjency[current_node]) > 0:
            current_node = adjency[current_node].pop(-1)
            continue
            
def FindPeptideRealSpectra(real_spectra):
    print("Make sure real spectra has node 0 added in the beginning")
    
    nodes = [i for i in range(len(real_spectra))]
    full_spectrum = nodes
    if 0 not in full_spectrum:
        full_spectrum.append(0)
    full_spectrum = sorted(full_spectrum)
    adj = Adjency(full_spectrum)
    start = nodes[0]
    end = nodes[-1]
    
    all_paths = FindAllPathWeightedDAG(start,end,adj,real_spectra)

    graph = Graph(full_spectrum)
    all_sums = []
    for path in all_paths:
        all_sums.append(sum([real_spectra[i] for i in path]))
    best_path = all_paths[np.argmax(all_sums)]

    peptide = ""
    for i in range(len(best_path)-1):
        prefix = full_spectrum[best_path[i]]
        suffix = full_spectrum[best_path[i+1]]
        edge = f"{prefix}->{suffix}"
        peptide += graph[edge]
    return peptide
    