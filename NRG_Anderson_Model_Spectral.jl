#################################################################################################################
#                                                 Importing                                                      #
#################################################################################################################

using LinearAlgebra  #Pacote de Algebra linear
using SpecialFunctions  #Pacote de Funções Especiais
using BitBasis, Base.Threads      # numbers in binary representation

#################################################################################################################
#                                                 Parâmetros                                                    #
#################################################################################################################

N_max=50            # Número de iterações
E_corte=20         # Energia de Corte
G=3.0              # Lambda
Ed=-15.0           # Energia do nível unicamente ocupado
U=30.0             # Energia de respulsão
GammaV=1.0         # Gamma   (proporcional ao quadrado da hibridização)


#################################################################################################################
#                                                 Function                                                      #
#################################################################################################################




#################################################################################################################
#                                          Função para Iteração                                                 #
#################################################################################################################

function Iteraction(Ed,U,GammaV,E_corte,G,N_max)
    
    file=open("EnergiasSIAM_Lambda3.txt","w")
    write(file, "###########################################################################################\n")
    write(file, "                E_d=",string(Ed),"  U=",string(U),"  Gama=",string(GammaV),"  Lambda=",string(G),"  E_corte=",string(E_corte),"\n")
    write(file, "###########################################################################################\n")
    write(file, "---------------------------------------------------------------------------\n---------------------------------------------------------------------------\n")  
        
    FacVivaldo=(G-1.0)/(G*log(G))      # Fator do Vivaldo

    GammaV=GammaV/G 
    Utio=U/FacVivaldo
    Edtio=Ed/FacVivaldo

    E_0=[0.0,Edtio,2.0*Edtio+Utio]

    E=(E_0.-minimum(E_0))./G
    println("Energias Iniciais ===> |vac>= ",E[1]," |+ ou ->= ",E[2]," |+,->= ",E[3],"\n")

    Regis=[-1 0; 0 +1; 1 0]            #Primeira coluna -> carga; Segunda coluna -> Spin

    # Elementos de matriz invariante

    ElemInv=zeros(3,3)
    ElemInv[2,1]=1.0
    ElemInv[3,2]=-sqrt(2.0)

    ###########################################################################################################################

    n=0
    ener=[]

    while n<=N_max

        write(file, "                                  N= ",string(n),"                                  \n")
               
        ContriVecN=Array{Int64}(undef, 0, 3)
        enerN_tot=[]
        ElemInvN=[]
        RegisN_tot=Array{Int64}(undef, 0, 2)
        Vec_eigN_tot=[]
        
        for q in (minimum(Regis[:,1])-1):1:(maximum(Regis[:,1])+1)
            for dS in (minimum(Regis[:,2])):1:(maximum(Regis[:,1])+1)
                    
                    m,s=NCadaGenero(Regis,q,dS)
                    
                    soma=sum(m)
                                    
                    if soma!=0
                                    
                        mH,VecContri=ConstMatHN(G,n,Regis,E,ElemInv,q,dS,GammaV)
                        
                        ContriVecN=vcat(ContriVecN, [q dS [VecContri]])
                        
                        E_eig=eigvals(mH)
                        Vec_eig=eigvecs(mH)
                                        
                        enerN_tot=vcat(enerN_tot,E_eig)
                        
                        for i in 1:1:size(E_eig,1)
                                
                                RegisN_tot=vcat(RegisN_tot, [q dS])   # adiciona uma ultima linha a RegisN com os valores de q e dS
                                push!(Vec_eigN_tot, Vec_eig[:,i])
                                
                        end
                                        
                    end
            end
        end
            
        E,Vec_eigN,RegisN=CutandPrint(E_corte,enerN_tot, Vec_eigN_tot, RegisN_tot,file)
        
        ElemInv=ElemInvNcalc(Vec_eigN,RegisN,Regis,ContriVecN)
        
        Regis=RegisN
            
        push!(ener,E)
            
        println("N: ",n," Dim: ",size(RegisN[:,1],1))
        
        write(file, "---------------------------------------------------------------------------\n---------------------------------------------------------------------------\n")
        write(file, "Dim= ",string(size(RegisN[:,1],1)),"\n")
        write(file, "---------------------------------------------------------------------------\n---------------------------------------------------------------------------\n")
        
        n+=1
               
        println("---------------------------------------------------------\n---------------------------------------------------------")
    end

    close(file)

    return "Fim"
end