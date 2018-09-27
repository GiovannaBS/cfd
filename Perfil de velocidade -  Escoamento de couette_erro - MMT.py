import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time

#Codigo modificado de escoamento entre duas placas estacionárias para escoamento de Couette
#OBSERVAÇÕES: Posteriormente, calcular P para qualquer dP e L

"""print("Perfil de velocidade do escoamento de couette")

print ("Densidade = 1000 (kg/m³) ")
print ("Viscosidade = 0.00101 (Ns/m²) ")
print ("Diâmetro = 0.02 (m) ")
print ("Coeficiente P = -1 ")
print('',end='\n')"""

ro = 1000
mi = 0.00101
D = 0.02
P = 1
E = 1

Re1 = 500
#Re2 = 1000
#Re3 = 2000

#t_test1 > t_test2
t_test1 = 100
t_test2 = 500

n_cels1 = 25
n_cels2 = 100
n_cels3 = 25

#print("n_cels1 ",n_cels1)
#print("n_cels2 ",n_cels2)
#print("n_cels3 ",n_cels3)

n_cels = n_cels1
Re=Re1

dy = 1/n_cels1
dt = E*Re*dy*dy
#dt=0.1
step_t = int(t_test2/dt)
print("step_t",step_t)

dy = 1/n_cels2
dt = E*Re*dy*dy
#dt=0.1
step_t1 = int(t_test2/dt)
print("step_t1",step_t1)

dy = 1/n_cels3
dt = E*Re*dy*dy
#dt=0.1
step_t2 = int(t_test2/dt)
print("step_t2",step_t2)

dy=0
dt=0


"""
for h in range (0,3,1):
 """       
t=t_test1
ini = time.time()
        
for d in range (0,2,1):
               
        for k in range (0,3,1):

                dy = 1/n_cels
                #print('dy',dy)

                dt = E*Re*dy*dy #definindo o intevalo dt - condição cfl para o metodo de crank-nicolson
                #print('dt',dt)
                #dt=0.1
#Utilizou-se Diferenças Finitas -  Método de Crank-Nicolson
#OBS.: Não entendi porque a constante 2P/Re não foi multiplicada por dt
#como eu perguntei antes

                s = dt/(2*Re*round(dy**2,6))
                        #print ("s",s)
                        #print ("2P/Re",(2*P/Re))
                        #print ("1-2s",1-2*s)
                        #print ("1+2s",1+2*s)
                        #print('',end='\n')

#Todas as variáveis foram adimensionadas

                u0 = 0 # velocidade u para y = 0 - condiçao de contorno (u0)
                u1 = 1 # velocidade u para y = d - condiçao de contorno (u1)
                t1 = int(t/dt) #numero de intervalos de tempo
                print("t1",t1)
                aux1=t1

# Matriz das velocidades (t1+1)x(n_cels+1)- definidas as condições de contorno
#o numero de linhas e colunas foi definido pelo numero d
#e intervalos do espaço (define o numero de colunas - incognitas das velocidades)
#e do tempo (define o numero de linhas - as novas velocidades no tempo n+1)

                u = np.zeros ((t1+1,n_cels+1))
                u[:,n_cels] = 1

#MÉTODO DA MATRIZ TRIDIAGONAL - para resolver o sistema
#Tenho n_cels intervalos, n_cels+1 elementos, dos quais 2 são condições de contorno
        
                a = np.zeros(n_cels-1)                  
                for i in np.arange(1,n_cels-1,1):
                        a[i] = -s
                b = np.zeros(n_cels-1)
                for i in np.arange(1,n_cels,1):
                        b[i-1] = 1+2*s
                c = np.zeros(n_cels-1) 
                for i in np.arange(1,n_cels-1,1):
                        c[i-1] = -s

#Calculo de cada linha para encontrar todas as incógnitas no tempo n+1
#utilizando a linha anterior de velocidades no tempo n
#começo calculando a linha j=1 (tempo t = 1/tempo anterior t = 0 -
#tenho a linha 0 com as primeiras velocidades para o tempo inicial)
#no python, pelo menos, se faz as contagens dos elementos de uma matriz começando de 0 a n-1
#nesse caso, estou me referindo a n como o numero de elementos da matriz - so como exemplo mesmo
#não é o n utilizado no código.

                for j in np.arange(1,t1+1,1):
    #o comando np.arange foi utilizado para gerar um array.
    #Antes tinha utilizado np.range, mas gerava uma lista e eu não podia usar no calculo com matriz
    #Explicando a contagem do comando - np.arange(>=limite inferior,<limite superior,passo)

                        r = np.zeros(n_cels-1) 
                        for i in np.arange(1,n_cels,1):
                                r[i-1] = s*(u[j-1,i+1] + u[j-1,i-1])+((1-2*s)*u[j-1,i])+(2*P*dt)/Re
           

    #Devido as condiçoes de contorno ja conhecidas, deixo em evidencia apenas as incognitas
    #Assim, realizo essa soma no primeiro e ultimo termo de cada intervalo
    #Ver Fortuna, a partir da pág 103 
                        r[0] = s*u0 + r[0]
                        r[n_cels-2] = s*u1 + r[n_cels-2]
        
#Resolvendo o sistema da matriz tridiagonal - link da referencia em anexo

                        if j==1:
                                c1 = np.zeros(n_cels-1)
                                c1[0]=c[0]/b[0]
                                for i in np.arange(1,n_cels-2,1): 
                                        c1[i]=c[i]/(b[i]-a[i]*c1[i-1])  
        
                        r[0]=r[0]/b[0]  
                        for i in np.arange(1,n_cels-1,1):   
                                r[i]=(r[i]-a[i]*r[i-1])/(b[i]-a[i]*c1[i-1])  
             
                        x=np.zeros(n_cels-1)
                        x[n_cels-2]=r[n_cels-2]
                        for i in np.arange(n_cels-3,-1,-1):
                                x[i] = r[i]-c1[i]*x[i+1]
                
                        for i in np.arange (1,n_cels,1):
                                u[j,i]=x[i-1]
                        
                #print ("nova matriz u")
                #print(u)
                #print(u.shape)
                #print('',end='\n')

#Ao fim, terei a matriz das velocidades atualizada com uma nova linha
#O processo se reenicia para o calculo da nova linha das velocidades no tempo n+1
#Utilizando sempre a linha anterior

                y = np.zeros(n_cels+1)                  
                for i in np.arange(1,n_cels+1,1):
                        y[i] = y[i-1]+dy
                #print ('y', y.shape)


                if k==0:
                        u_test16=u
                        y_test16=y
                        n_cels=n_cels2
                if k==1:
                        u_test17=u
                        y_test17=y
                        n_cels=n_cels3
                if k==2:
                        u_test18=u
                        y_test18=y
                        n_cels=n_cels1
                        
        if d==0:
                t=t_test2
                u_test13=u_test16
                y_test13=y_test16
                u_test14=u_test17
                y_test14=y_test17
                u_test15=u_test18
                y_test15=y_test18
                fim = time.time()
                ini1 = time.time()
fim1 = time.time()
print ("Tempo de execução para "+str(t_test1)+": ", fim-ini)
print ("Tempo de execução para "+str(t_test2)+": ", fim1-ini1)
                        
"""
if h==0:
        Re=Re2"""
        #para Re=Re1/t=t1
u_test1 = u_test13
y_test1=y_test13
u_test2 = u_test14
y_test2=y_test14
u_test3 = u_test15
y_test3=y_test15
        #para Re=Re1/t=t2
u_test4 = u_test16
y_test4=y_test16
u_test5 = u_test17
y_test5=y_test17
u_test6 = u_test18
y_test6=y_test18

""" 
if h==1:
        Re=Re3
        #para Re=Re2/t=t1
        u_test7 = u_test13
        y_test7=y_test13
        u_test8 = u_test14
        y_test8=y_test14
        u_test9 = u_test15
        y_test9=y_test15
        #para Re=Re2/t=t2
        u_test10 = u_test16
        y_test10=y_test16
        u_test11 = u_test17
        y_test11=y_test17
        u_test12 = u_test18
        y_test12=y_test18
"""          
                        
#SOLUÇÃO ANALITICA

u_an0 = np.zeros((1,n_cels+1))
y_an0 = np.zeros(n_cels+1)
for i in np.arange(0,n_cels+1,1):  
        if i>0:
                y_an0[i] = y_an0[i-1]+(1/n_cels)
        u_an0[0,i] = y_an0[i]+(P*y_an0[i]*(1-y_an0[i]))



#ERROS
#print ("n_cels1", n_cels1)

#erro para tds as n células
erro1 = np.zeros ((step_t+1,n_cels1+1))
erro2 = np.zeros ((step_t1+1,n_cels2+1))
erro3 = np.zeros ((step_t2+1,n_cels3+1))

#erro máximo observado entre tds as células em cada passo de tempo
e_max1 = np.zeros ((step_t+1,1))
e_max2 = np.zeros ((step_t1+1,1))
e_max3 = np.zeros ((step_t2+1,1))

"""arrays para solução analitica foram criados porque as matrizes de velocidade 
tem n colunas - variavel função do numero de células - a fim de que a solução
analitica possua mesmo numero de colunas
o objetivo é gerar 3 gráficos (um para cada Re) e observar a variação do erro
maximo no tempo em função também do numero de células"""

u_an1 = np.zeros((1,n_cels2+1))
y_an1 = np.zeros(n_cels2+1)
for i in np.arange(0,n_cels2+1,1):  
        if i>0:
                y_an1[i] = y_an1[i-1]+(1/n_cels2)
        u_an1[0,i] = y_an1[i]+(P*y_an1[i]*(1-y_an1[i]))

u_an2 = np.zeros((1,n_cels3+1))
y_an2 = np.zeros(n_cels3+1)
for i in np.arange(0,n_cels3+1,1):  
        if i>0:
                y_an2[i] = y_an2[i-1]+(1/n_cels3)
        u_an2[0,i] = y_an2[i]+(P*y_an2[i]*(1-y_an2[i]))

"""print ('u_an0', u_an0.shape)
print ('u_an1', u_an1.shape)
print ('u_an2', u_an2.shape)
print ('erro1', erro1.shape)
print ('erro2', erro2.shape)
print ('erro3', erro3.shape)
print ('u_test4', u_test4.shape)
print ('u_test5', u_test5.shape)
print ('u_test6', u_test6.shape)
"""
for j in np.arange(0,step_t+1,1):
        for i in np.arange(0,n_cels1+1,1):       
                #Re = 500
                erro1[j,i]=abs(u_test4[j,i]-u_an0[0,i])
                if i == 0:
                        e_max1[j,0] = erro1[j,i]
                else:
                        if e_max1[j,0] < erro1[j,i]:
                                e_max1[j,0] = erro1[j,i]
        if j>step_t-3:
                print ('e_max1',e_max1[j,0])
                                
                #if i <= n_cels2:
                    #Re = 500
                        #erro2[j,i]=abs(u_test5[j,i]-u_an1[0,i])
                        #if i == 0:
                                #e_max2[j,0] = erro2[j,i]
                        #else:
                                #if e_max2[j,0] < erro2[j,i]:
                                        #e_max2[j,0] = erro2[j,i]
                #if i <= n_cels3:
                #Re = 500
                        #erro3[j,i]=abs(u_test6[j,i]-u_an2[0,i])
                        #if i == 0:
                                #e_max3[j,0] = erro3[j,i]
                        #else:
                                #if e_max3[j,0] < erro3[j,i]:
                                        #e_max3[j,0] = erro3[j,i]

for j in np.arange(0,step_t1+1,1):
        for i in np.arange(0,n_cels2+1,1): 
                #Re = 500
                erro2[j,i]=abs(u_test5[j,i]-u_an1[0,i])
                if i == 0:
                        e_max2[j,0] = erro2[j,i]
                else:
                        if e_max2[j,0] < erro2[j,i]:
                                e_max2[j,0] = erro2[j,i]
        if j>step_t1-3:
                print ('e_max2',e_max2[j,0])
                        
for j in np.arange(0,step_t2+1,1):
        for i in np.arange(0,n_cels3+1,1):   
                #Re = 500
                erro3[j,i]=abs(u_test6[j,i]-u_an2[0,i])
                if i == 0:
                        e_max3[j,0] = erro3[j,i]
                else:
                        if e_max3[j,0] < erro3[j,i]:
                                e_max3[j,0] = erro3[j,i]
        if j>step_t2-3:
                print ('e_max3',e_max3[j,0])

del erro1, erro2, erro3  
                                              
"""
plt.figure()
matplotlib.rcParams.update({'font.size': 6})

#ou font = {'family' : 'normal','weight' : 'bold','size'   : 22}
#   matplotlib.rc('font', **font)

plt.subplot(1,1,1)
plt.ylabel('y/D', fontsize=10)
plt.xlabel('u/U', fontsize=10)
plt.grid(linestyle = ':')
plt.title('Re = 2000/ t = 100/ P = -1', fontsize=11)
plt.plot(u_test1[-1,:],y_test1,'m',label = "n_cels = "+str(n_cels1)+"")  #-1 mostra a ultima linha da matriz
plt.plot(u_test2[-1,:],y_test2,'c',label = "n_cels = "+str(n_cels2)+"")
plt.plot(u_test3[-1,:],y_test3,'r',label = "n_cels = "+str(n_cels3)+"")
plt.plot(u_an0[0,:],y_an0,'k',label = "Solução analítica")
plt.tick_params(axis='both', labelsize=10)
plt.legend(loc = 'upper left', prop={'size':10})
#plt.savefig('Re_2000_t_100_pm1.png')
#plt.show()


plt.subplot(1,1,1)
plt.ylabel('y/D', fontsize=10)
plt.xlabel('u/U', fontsize=10)
plt.grid(linestyle = ':')
plt.subplots_adjust(hspace=0.5)
#subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
plt.title('Re = 2000/ t = 500/ P = -1', fontsize=11)
plt.plot(u_test4[-1,:],y_test4,'m',label = "n_cels = "+str(n_cels1)+"")  #-1 mostra a ultima linha da matriz
plt.plot(u_test5[-1,:],y_test5,'c',label = "n_cels = "+str(n_cels2)+"")
plt.plot(u_test6[-1,:],y_test6,'r',label = "n_cels = "+str(n_cels3)+"")
plt.plot(u_an0[0,:],y_an0,'k',label = "Solução analítica")
plt.tick_params(axis='both', labelsize=10)
plt.legend(loc = 'upper left', prop={'size':10})
#plt.savefig('Re_2000_t_500_pm1.png')
#plt.show()


step =np.zeros(step_t+1)
for i in np.arange(0,step_t+1,1):  
        step[i] = i/step_t
step1 =np.zeros(step_t1+1)
for i in np.arange(0,step_t1+1,1):  
        step1[i] = i/step_t1
step2 =np.zeros(step_t2+1)
for i in np.arange(0,step_t2+1,1):  
        step2[i] = i/step_t2


plt.subplot(1,1,1)
plt.xlabel('i/step', fontsize=10)
plt.ylabel('Erro máximo', fontsize=10)
plt.grid(True, linestyle = ':',which="both")
plt.title('Erro máximo (Re = 500/ t = 500/ P = 1)', fontsize=11)
plt.plot(step,e_max1[:,0],'m',label = "n_cels = "+str(n_cels1)+"")  #-1 mostra a ultima linha da matriz
plt.plot(step1,e_max2[:,0],'c',label = "n_cels = "+str(n_cels2)+"")
plt.plot(step2,e_max3[:,0],'r',label = "n_cels = "+str(n_cels3)+"")
plt.xscale('log')
plt.yscale('log')
plt.tick_params(axis='both', labelsize=10)
plt.legend(loc = 'lower left', prop={'size':10})
#plt.savefig('Re_500_t_500_erro_p1.png')
#plt.show()

"""

