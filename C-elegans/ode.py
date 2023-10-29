# Created by Kangxin Hu at Beihang University, on Oct.2023.

def odesolve(the_list,y,y0):
    y=[0 for i in range(0,len(the_list))]
    y[0]=y0
    for j in range(1,len(the_list)):
        y[j]=y[j-1]+0.02*(-the_list[j-1]*0.001-y[j-1]*0.2)/0.1
    return y

    with open("data1.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[278]))
    with open("data2.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[279]))
    with open("data3.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[280]))
    with open("data4.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[281]))
    with open("data5.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[282]))
    with open("data6.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[283]))
    with open("data7.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[284]))
    with open("data8.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[285]))
    with open("data9.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[286]))
    with open("data10.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[287]))
    with open("data11.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[288]))
    with open("data12.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[289]))
    with open("data13.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[290]))
    with open("data14.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[291]))
    with open("data15.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[292]))
    with open("data16.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[293]))
    with open("data17.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[294]))
    with open("data18.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[295]))
    with open("data19.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[296]))
    with open("data20.csv", "w") as f:
        csv.writer(f).writerows(zip(t, neurons_v[297]))




