#temperature of cool surface(first assume)
t_sc=53.5


##linear interpolation a=value, p=properties dictionary
def linear(a,p):
  key1=a//5*5
  key2=(a//5+1)*5
  result=(p[key1]+((p[key2]-p[key1])/5)*(a-key1))
  return result

##properties dictionary
density={15:999.1,20:998,25:997,30:996,35:994,40:992.1,45:990.1,
     50:988.1,55:985.2}
viscosity={15:1.138*10**-3,20:1.002*10**-3,25:0.891*10**-3,30:0.798*10**-3
           ,35:0.72*10**-3,40:0.653*10**-3,45:0.596*10**-3,50:0.547*10**-3
           ,55:0.504*10**-3}
pr={15:8.09,20:7.01,25:6.14,30:5.42,35:4.83,40:4.32,45:3.91,50:3.55}
beta={15:0.138*10**-3,20:0.195*10**-3,25:0.247*10**-3,30:0.294*10**-3,
      35:0.337*10**-3,40:0.377*10**-3,45:0.415*10**-3,50:0.451*10**-3}
kf={15:0.589,20:0.598,25:0.607,30:0.615,35:0.623,40:0.631,45:0.637,50:0.644}

error=1
while error>0.01 or error<-0.01:
  #film temp.
  t_f=(t_sc+7)/2
  #properties @ Tf
  r=linear(t_f,density)
  y=linear(t_f,viscosity)
  p=linear(t_f,pr)
  b=linear(t_f,beta)
  k=linear(t_f,kf)
  v=y/r
  gr=(9.81*b*(t_sc-7)*(0.2**3))/((v)**2)
  ra=gr*p
  nu=(0.825+(0.387*(ra**(1/6)))/((1+(0.492/p)**(9/16))**(8/27)))**2
  h=(nu*k)/0.2
  t_sc_new=((1500/0.025+7*h)/(h+15/0.025))
  error=((t_sc_new-t_sc)/t_sc_new)
  t_sc=t_sc_new
  print(t_sc)
  print(error)
