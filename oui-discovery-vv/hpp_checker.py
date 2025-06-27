from explore_import import  *

def get_pep_positions(prot,pep):    
    return [(m.start(),m.start()+len(pep)-1) for m in re.finditer(pep, prot)]
def check_ov_any(pep_pos,all_pep_pos):
    """ which of the positions of pep overlap with any other pep
    pep_pos is a list
    """
    for pos in pep_pos: all_pep_pos.remove(pos)
    other_pep_pos=sorted(all_pep_pos,key=lambda x: x[1])
    pep_pos_ov=[]
    for j,pos in enumerate(pep_pos):
        #print(pos)
        x=interval.interval[*pos]
        not_ov=True; i=0
        while not_ov and i<len(other_pep_pos): 
            y=interval.interval[*other_pep_pos[i]]
            if len(x & y):
                not_ov=False
            #print(not_ov)
            i+=1
        if not not_ov: pep_pos_ov=pep_pos_ov+[j]
    return pep_pos_ov
    
def check_type_len_ov(pep_pos,all_pep_pos,typeof_ov):
    """pep_pos is a tuple
        all_pep_pos is a list of tuples, not sorted
        typeof_ov="partial" or "total"
    """
    x=interval.interval[pep_pos]
    typeof_ov_ifpass=[]; lenof_ov_ifpass=[]
    for other_pos in all_pep_pos:
        if other_pos==pep_pos: continue
        y=interval.interval[other_pos]
        ##if len(x & y)==0: return False, None
        if typeof_ov=="total":
            if len(x & y)>0 and  ((x & y)==x or (x & y)==y): #two-sided                
                max_len=np.max([pep_pos[1]-pep_pos[0],other_pos[1]-other_pos[0]])
                res=(True, max_len>=9)
            else:
                res=(False, False) #None
        if typeof_ov=="partial":
            #print(x,y)
            if  len(x & y)>0 and ((x & y)!=x or (x & y)!=y): #partial overlap will yield false
                ov_int=sorted([pep_pos,other_pos],key=lambda x: x[1])
                tot_len=ov_int[1][1]-ov_int[0][0]
                #print(0)
                res=(True, tot_len>=18)
            else: 
                #print(1)
                res=(False, False) #None
        typeof_ov_ifpass=typeof_ov_ifpass+[res[0]]
        lenof_ov_ifpass=lenof_ov_ifpass+[res[1]]
    return  typeof_ov_ifpass, lenof_ov_ifpass
        
def HPP_decision(pep,pep_pos_ov,type_len_ov_any):
    """returns if info has to be collected, and if msms counts should be merged to other pep"""
    if len(pep_pos_ov)==0:
        return len(pep)>=9, False
    else:
        if type_len_ov_any==[True,False,True,True] or type_len_ov_any==[False,None,True,True]:
            return True, False
        else:
            #if len(pep)>=9:
            #    return True, False
            #else:
            if type_len_ov_any[1]==True or type_len_ov_any[3]==True:
                return True, True #add msms count to 1 or both of overlap peps
            else:
                return False, False#"mess with code" #is it???????????????????????                 
def find_max_length_in_overlap_groups(intervals):
    #https://stackoverflow.com/questions/77116236/find-longest-interval-between-overlapping-intervals-in-python
    #store the indices
    indices={val:ind for ind, val in enumerate(intervals)}
    aux = sorted(intervals, key=lambda x: x[0])
    #print(aux)
    new_group = [aux[0]]
    output = []
    for interval in aux[1:]:
        if interval[0] > new_group[-1][-1]:
            output.append(new_group)
            new_group = [interval]
        else:
            new_group.append(interval)
    output.append(new_group)
    #print(output)
    elementindx_intervals=[[indices[i] for i in intervals] for intervals in output]
    longest_intervals=[(np.min([interval[0] for interval in intervals]),np.max([interval[1] for interval in intervals])) for intervals in output]
    #print(longest_intervals)
    #print(elementindx_intervals)
    #sort the output according to indices in input to preserve order
    return elementindx_intervals, longest_intervals