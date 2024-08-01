#################################################################################################################
#                                                 Parâmetros                                                    #
#################################################################################################################

N_max=5            # Número de iterações
E_corte=5         # Energia de Corte
G=1.0              # Lambda
Ed=-0.0           # Energia do nível unicamente ocupado
U=0.0             # Energia de respulsão
V=0.1         # Gamma   (proporcional ao quadrado da hibridização)

#################################################################################################################
#                                                 Importing                                                      #
#################################################################################################################

using LinearAlgebra  #Pacote de Algebra linear
using SpecialFunctions  #Pacote de Funções Especiais
using BitBasis, Base.Threads      # numbers in binary representation

#################################################################################################################
#                                                 Function                                                      #
#################################################################################################################

#################################################################################################################
# Determina o numero de estados de um certo genero, g, dados: o Registro (Regis) da iteração N-1, a carga e o spin 
#################################################################################################################

function NGenero(Regis,g,Q,dS)
    
    numb=0
    posi=0
    lug=[]
    
    if g==1  #Norte
        
        for i in 1:1:size(Regis,1)
                       
            posi+=1
            
            if (Regis[i,1]+1)==Q && Regis[i,2]==dS
                numb+=1
                push!(lug,posi)
            end
            
        end
        
    elseif g==2  #Sul
        
        for i in 1:1:size(Regis,1)
                       
            posi+=1
            
            if (Regis[i,1]-1)==Q && Regis[i,2]==dS
                numb+=1
                push!(lug,posi)
            end
            
        end
        
    elseif g==3    # Leste

        for i in 1:1:size(Regis,1)

            posi+=1
            
            if Regis[i,1]==Q && (Regis[i,2]+1)==dS
                numb+=1
                push!(lug,posi)
            end
            
        end
        
    elseif g==4    # Oeste
                
        for i in 1:1:size(Regis,1)

            posi+=1
            
            if Regis[i,1]==Q && (Regis[i,2]-1)==dS
                numb+=1
                push!(lug,posi)
            end
            
        end
        
    end
        
    
    return numb, lug
end

#################################################################################################################
# Determina o numero de contribuições de N, S, L e W para cada iteração e os estados que geram tais gêneros
#################################################################################################################

function NCadaGenero(Regis,Q,dS)
    
    mat=[]  # numero de estados com certo gênero
    estat=[]   # estados que contribuem em cada gênero
    
    # 1 -> Norte
    # 2 -> Sul
    # 3 -> Leste
    # 4 -> Oeste
    
    for i in 1:4
        
        numb,lug=NGenero(Regis,i,Q,dS)
        
        push!(mat,numb)
        
        if lug!=[]
            
            for i in lug
                push!(estat,i)
            end
            
        end
            
    end
    
    return mat, estat
end

#################################################################################################################
# Controi a matrix do Hamiltoniano na iteração N por sertor de carga e spin
#################################################################################################################

function ConstMatHN(G,N,Regis,E,ElemInv,Q,dS,V)
    
    GeneroContri,VecContri=NCadaGenero(Regis,Q,dS)   # Encontra as direções (gêneros) contribuintes e os autoestados contribuintes
    
    dim=sum(GeneroContri)      # encontra a dimensão deste setor de H de carga Q e spin dS
    
    mH=zeros(dim,dim)          # inicializa a matriz
    
    # Elementos fora da diagonal
    
    for i in 1:1:dim
        for j in 1:1:dim
            
            mH[i,j]+=ElemMatH(i,j,G,N,dS,Regis,ElemInv,GeneroContri,VecContri,V)
            
        end
    end

    mH1=mH+transpose(conj(mH)) # Tornando hermitiano
     
    for i in 1:1:dim
        
        mH1[i,i]+=sqrt(G)*E[VecContri[i]]
        
    end
       
    return mH1,VecContri
end

#################################################################################################################
# Acha o genero de um dado elemento da base de H_N
#################################################################################################################

function FindGenN(i,GeneroContri)
    
    if i<=GeneroContri[1]
        
        return 1     # Norte
        
    elseif GeneroContri[1]<i<=(GeneroContri[2]+GeneroContri[1])
        
        return 2     # Sul 
        
    elseif (GeneroContri[2]+GeneroContri[1])<i<=(GeneroContri[2]+GeneroContri[1]+GeneroContri[3])
        
        return 3     # Leste
        
    elseif (GeneroContri[2]+GeneroContri[1]+GeneroContri[3])<i<=(GeneroContri[2]+GeneroContri[1]+GeneroContri[3]+GeneroContri[4])
    
        return 4     # Oeste
        
    end

end

#################################################################################################################
# Vê se os estados são compatíveis pelo gênero na iteração N
#################################################################################################################

function CompOrNotComp(i,j,GeneroContri)
    
    mat=[0 0 1 1; 0 0 1 1; 1 1 0 0; 1 1 0 0] # matriz de não nulidade dos elementos
    
    if (mat[FindGenN(i,GeneroContri),FindGenN(j,GeneroContri)]==0)
        
        return "NotComp" # Não compatível
        
    elseif mat[FindGenN(i,GeneroContri),FindGenN(j,GeneroContri)]==1

        return "Comp"   # Compatível
        
    end
        
end

#################################################################################################################
#Calcula E_{N-1}
#################################################################################################################

 function ENless1(G,N)
    
    en=(1-G^(-1.0))/log(G)*G^(-1.0*N)
    
    return en
end

#################################################################################################################
# Calcula os elementos de matriz do Hamiltoniano
#################################################################################################################

function ElemMatH(i,j,G,N,dS,Regis,ElemInv,GeneroContri,VecContri,GammaV)
    
    el=0.0
    
    if CompOrNotComp(i,j,GeneroContri)=="NotComp"  # Os estados são não compativeis
        
        el+=0.0    # logo geram elemento nulo
        
    else      # se são compativeis, calcule o elemento
        
          el+=tENless1(G,N,GammaV)*Braf_dagger_fKet(i,j,dS,ElemInv,GeneroContri,VecContri)
                    
    end
        
    return el
    
end

#################################################################################################################
# Calcula o elemento de matriz de f^{dagger}f
#################################################################################################################

function Braf_dagger_fKet(i,j,dS,ElemInv,GeneroContri,VecContri)
    
    mat=[0.0 0.0 sqrt(dS/(dS+1.0)) -sqrt((dS+2.0)/(dS+1.0)); 0.0 0.0 1.0 1.0; sqrt(dS/(dS+1.0)) 1.0 0.0 0.0; -sqrt((dS+2.0)/(dS+1.0)) 1.0 0.0 0.0]
    
    # Gênero dos estados i e j
    
    gi=FindGenN(i,GeneroContri)
    gj=FindGenN(j,GeneroContri)
    
    # Pais dos estados i e j
    
    Pi=Int(VecContri[i])
    Pj=Int(VecContri[j])
    
    return  mat[gi,gj]*ElemInv[Pi,Pj]
    
end

#################################################################################################################
# Calculo dos elementos de matriz invariante
#################################################################################################################

function ElemInvNcalc(Vec,RegisN,Regis,ContriVecN)
    
    dim=size(RegisN[:,1],1)
        
    matElemInvN=zeros(dim,dim)
  
   @threads for i in 1:1:dim
        for j in 1:1:dim

            if i!=j
            
                matElemInvN[i,j]+=InvVecs(i,j,Vec,RegisN,ContriVecN,Regis)

            end
            
        end
    end
    
    return matElemInvN
end

function InvVecs(i,j,Vec,RegisN,ContriVecN,Regis)
    
    Qi=RegisN[i,1]
    dSi=RegisN[i,2]
    Qj=RegisN[j,1]
    dSj=RegisN[j,2]
    
    el=0.0
       
    if (Qi-Qj)==1 && abs(dSi-dSj)==1
             
        ContriVeci=ReadQdSContri(Qi,dSi,ContriVecN)
        ContriVecj=ReadQdSContri(Qj,dSj,ContriVecN)
        
        Veci=Vec[i]
        Vecj=Vec[j]
        
        for p in 1:1:size(Veci,1)
            
            Pp=ContriVeci[p]
            Gp=FindGPN(Qi,dSi,p,Regis)
            
            for pl in 1:1:size(Vecj,1)
                
                Ppl=ContriVecj[pl]
                Gpl=FindGPN(Qj,dSj,pl,Regis)
                Mgg=Mggl(dSj,Int(Gp),Int(Gpl))
                
                if Mgg!=0.0 && Ppl==Pp
                    
                    el+=Veci[p]*Mgg*Vecj[pl]
                    
                end
            end
            
        end
    end
    
    return el
    
end

###########################################################################################################################
# Encontra o gênero e o estado pai de uma certa posição de um dado autovetor do setor de carga Q e spin dS do Hamiltoniano
###########################################################################################################################

function FindGPN(Q,dS,i,Regis)   # Regis é o matriz regitro da iteração N-1
    
    NcadaGen, EstatPai=NCadaGenero(Regis,Q,dS)
    
    Gen_i=FindGenN(i,NcadaGen)
     
    return Gen_i
    
end

#################################################################################################################
# Retorna a Conexão M_{g,g'} dos estados com gêneros g e g'
#################################################################################################################

function Mggl(dS,g,gl)
   
    mat2=zeros(4,4)
    
    mat2[1,3]=-sqrt((dS+1)/dS)
    mat2[1,4]=sqrt((dS+1)/(dS+2))
    mat2[3,2]=1.0
    mat2[4,2]=1.0
    
                
    return mat2[g,gl]
        
end

#################################################################################################################
# Encontra, dado Q e dS, a posição inicial e o numero de posições do setor de carga Q e spin dS
#################################################################################################################

function InFimSecQdS(Q,dS,RegisN)   # RegisN é a matriz registro da iteração N
    
    inicio=0
    matcontrol=[]
    
    for i in 1:1:size(RegisN[:,1],1)
               
        inicio+=1
            
        if RegisN[i,1]==Q && RegisN[i,2]==dS
            
            push!(matcontrol,inicio)
            
        elseif  RegisN[i,1]!=Q && RegisN[i,2]!=dS && matcontrol!=[]    
            break
        end
        
    end
    
    if matcontrol!=[]
        return  Int(minimum(matcontrol)), Int(maximum(matcontrol))  # retorna o inicio do setor e o fim do setor, respectivamente
    else
        return 1,1
    end
    
end

#################################################################################################################
# Cortando e escrevendo
#################################################################################################################

function CutandPrint(E_corte,E_eig_tot, Vec_eig_tot, RegisN_tot,file)
    
    RegisN=Array{Int64}(undef,0, 2)
    enerN=[]
    Vec_eigN=[]
    
    enfund=minimum(E_eig_tot)
    
    println("Energia Fundamental: ", enfund)
    
    write(file, "---------------------------------------------------------------------------\n","Energia Fundamental: ", string(enfund),"\n---------------------------------------------------------------------------\n")
    write(file," Q   |   dS   ======>   E \n")

    for i in 1:1:size(RegisN_tot,1)
                                                   
        en=E_eig_tot[i]-enfund
                        
        if en<=E_corte+enfund            # Corta as energias e atualiza a lista de energia e o Registro
            
            push!(enerN,en)
            RegisN=vcat(RegisN,[RegisN_tot[i,1] RegisN_tot[i,2]])   # adiciona uma ultima linha a RegisN com os valores de q e dS
            push!(Vec_eigN, Vec_eig_tot[i])
            #println(RegisN_tot[i,1],' ',RegisN_tot[i,2]," =====> ",en)
            
            write(file," ",string(RegisN_tot[i,1])," | ",string(RegisN_tot[i,2])," ======> ",string(en),"\n")
    
        end
    end
    
    println("---------------------------------------------------------")
    
    return enerN, Vec_eigN, RegisN
    
end               

#################################################################################################################
# Le as contrituições para a carga Q e spin dS
#################################################################################################################

function ReadQdSContri(Q,dS,ContriVecN)
    
    g=[]
    
    for i in 1:1:size(ContriVecN,1)
        
        q=ContriVecN[i,1]
        ds=ContriVecN[i,2]
        
        if  q==Q && ds==dS
            g=ContriVecN[i,3]
        end
    end
    
    return g
    
end

#################################################################################################################
# Calcula t_{N-1}/E_{N-1}
#################################################################################################################

function tENless1(G,N,GammaV)
    
    tn=0.0
    #FacVivaldo=(G-1.0)/(G*log(G))
    
    if N==0
                
        #GammaVtio=2.0*GammaV/pi
        
        tn+=sqrt(2)*V
        
    else
        N-=1
       tn+=G^(-N-1/2)
         
    end
    
    return tn
    
end

#################################################################################################################
# D_N                                                                                                           #
#################################################################################################################

function D_N(G,N)
    
    return (G-1)/(G*log(G))*G^(-1.0*(N-1.0)/2)
    
end

#################################################################################################################
#                                          Função para Iteração                                                 #
#################################################################################################################

function Iteraction(Ed,U,V,E_corte,G,N_max)
    
    file=open("eNRG_lamb1_EnergiasSIAM.txt","w")
    write(file, "###########################################################################################\n")
    write(file, "                E_d=",string(Ed),"  U=",string(U),"  Gama=",string(V),"  Lambda=",string(G),"  E_corte=",string(E_corte),"\n")
    write(file, "###########################################################################################\n")
    write(file, "---------------------------------------------------------------------------\n---------------------------------------------------------------------------\n")  
        
    #FacVivaldo=(G-1.0)/(G*log(G))      # Fator do Vivaldo

    V=V
    Utio=U
    Edtio=Ed

    Euns=[1.0,1.0,1.0]
    E_0=[0.0,Edtio,2.0*Edtio+Utio]

    E=(E_0-Euns*minimum(E_0))/G
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
                                    
                        mH,VecContri=ConstMatHN(G,n,Regis,E,ElemInv,q,dS,V)
                        
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

#################################################################################################################
#                                                 Rodando                                                       #
#################################################################################################################

@time Iteraction(Ed,U,V,E_corte,G,N_max)