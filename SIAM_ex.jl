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

V=0.1                               # Hybridization
U=2.0                               # Repulsion
E_d=U/2                            # Impurity Energy
L=5                                 # Bumber of Bath sites

############################################################################################################################
#################################################---- Functions ----########################################################
############################################################################################################################


#######################################---- Generating Basis States ----#####################################################

function Basis(L)

    local_states=["vac","up","down","updown"]

    States=[["vac"],["up"],["down"],["updown"]]



    for l=2:1:L

        sup=[]

        for st in States

            vec=copy(st)

            for j=1:4

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

    if local_state=="vac"

        return 0.0

    elseif local_state=="up" || local_state=="down"

        return E_d

    else

        return 2.0*E_d+U

    end

end

function  off_diagonal_mH(State1,State2,V)


    dif=findall(in(State1),State2) # Find where the vectors are iqual

    if abs(length(State1)-length(dif))==2

        if State1[1]==State2[1] && State1[2]==State2[2]

            return -1.0

        else

            return sqrt(2)*V

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


############################################################################################################################
#########################################---- Main Calculation function ----################################################
############################################################################################################################

function SIAM_ex(L,E_d,U,V)

    Basis_states=Basis(L)

    H=mH(Basis_states,E_d,U,V)

    E0=eigmin(H)

    return E0, H

end



############################################################################################################################
#################################################---- Really Calculation ----###############################################
############################################################################################################################

@time e0,h=SIAM_ex(L,E_d,U,V)