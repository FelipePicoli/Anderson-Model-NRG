############################################################################################################################
#################################################--- Photoemission ---######################################################
############################################################################################################################


############################################################################################################################
####################################################--- Importing ---#######################################################
############################################################################################################################

using LinearAlgebra  #Pacote de Algebra linear
using Random, Distributions  # random number's distributions
using BitBasis, Base.Threads      # numbers in binary representation
using LaTeXStrings  # For LaTex code

############################################################################################################################
#############################################--- Setting the parameters ---################################################
############################################################################################################################

V=0.1                              # Hybridization
U=5.0                               # Repulsion
E_d=-U/2                            # Impurity Energy
L=5                                 # Bumber of Bath sites
Occ_number=false                     # Do you want calculate the Occupation number?

############################################################################################################################
#################################################---- Functions ----########################################################
############################################################################################################################


#######################################---- Generating Basis States ----#####################################################

function Basis(L)

    #local_states=["vac","up","down","updown"]

    #States=[["vac"],["up"],["down"],["updown"]]

    local_states=[0,1,2]

    States=[[0],[1],[2]]


    for l=2:1:L

        sup=[]

        for st in States

            vec=copy(st)

            for j=1:1:length(local_states)

                push!(sup,push!(vec,local_states[j]))
                vec=copy(st)

            end

        end

        States=copy(sup)

    end

    return States

end

#######################################---- Calculating the Hamiltonian Elements ----###########################################

function  Diagonal_mH(State,E_d,U)

    local_state=State[1]

    if local_state==0

        return 0.0

    elseif local_state==1 

        return E_d

    elseif local_state==2

        return 2.0*E_d+U

    end

end

function  off_diagonal_mH(State1,State2,V)


    dif=findall(State1.!=State2) # Find where the vectors are iqual

    if length(dif)==2

        if abs(dif[1]-dif[2])==1 

            if ((State1[dif[1]]+State1[dif[2]])==(State2[dif[1]]+State2[dif[2]])) && abs(State1[dif[1]]-State2[dif[1]])==1 && abs(State1[dif[2]]-State2[dif[2]])==1

                if dif[1]==1 || dif[2]==1

                    return sqrt(2)*V
                else
                    return -1.0

                end

            else
                
                return 0.0
            end

        else
            return 0.0
        end

          
    else

        return 0.0

    end   
    
end

#############################################---- Calculating the Hamiltonian ----############################################

function mH(Basis_states,E_d,U,V)

    mat=zeros(length(Basis_states),length(Basis_states))

    for i=1:1:length(Basis_states)
        for j=i:1:length(Basis_states)

           if i==j
            
            mat[i,i]=Diagonal_mH(Basis_states[i],E_d,U)

           else

            mat[i,j]=off_diagonal_mH(Basis_states[i],Basis_states[j],V)
            mat[j,i]=mat[i,j]

           end

        end
    end

    return mat
end

###########################################---- Impurity Occupation Number ----############################################

function n_f(Vec,Basis_states)

    n=0.0

    for i=1:1:length(Vec)

        state=Basis_states[i]

        if state[1]==0

            N_el=0.0

        elseif state[1]==1 

           # N_el=1.0

            N_el=2.0

        elseif state[1]==2

            N_el=2.0

        end

        n+=abs2(Vec[i])*N_el

    end

    return n
end

############################################################################################################################
#########################################---- Main Calculation function ----################################################
############################################################################################################################

function SIAM_ex(L,E_d,U,V,Occ_number)

    Basis_states=Basis(L)

    H=mH(Basis_states,E_d,U,V)

    if Occ_number

        E,Vec=eigen(H)

        N_imp=n_f(Vec[:,1],Basis_states)

        return E[1],N_imp, H

    else

        E=eigmin(H)

        return E, H
    end

end



############################################################################################################################
#################################################---- Really Calculation ----###############################################
############################################################################################################################

if Occ_number
    @time e0,nf,H=SIAM_ex(L,E_d,U,V,Occ_number)
else
    @time e0,H=SIAM_ex(L,E_d,U,V,Occ_number)
end