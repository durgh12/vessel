import numpy as np
import math


# 주어진 상수 및 매개변수
E = 206000  # N/mm^2
v = 0.3
HY80orm = 552  # MPa
HY100orm = 686  # MPa
HY80 = HY80orm * 0.98
HY100 = HY100orm * 0.98
pi = np.pi
mildo = 7.85  # 단위 : ton/m**3
e = 2.7182818284
P = 5.5       # 설계압력
Rm = 3100     # 내경


# 설계변수 범위
h_range = np.arange(21,24, 0.5)             #외판두께
Lsp_range = range(370, 450, 10)             #늑골간격
hst_range = np.arange(160, 200, 5)          #웨브높이
tst_range = np.arange(18.0, 21.0, 0.5)      #웨브두께
bg_range = range(80, 100, 5)                #플랜지폭
tg_range = np.arange(31, 37, 0.5)           #플랜지 두께

# 결과를 저장 리스트
results = []

# 매개변수 조합을 위한 중첩된 루프
for h in h_range:
    for Lsp in Lsp_range:
        for hst in hst_range:
            for tst in tst_range:
                for bg in bg_range:
                    for tg in tg_range:
                        # 설계에 따른 bbtz 및 dbtz를 계산하는 수식 추가
                       
                        pii= 1.23* np.sqrt(Rm*h)/(Lsp -tst)

                        ohmx= -Rm*P/(2*h)
                        Asp= bg*tg+tst*hst
                        Gy=(bg*tg*(tg/2 +hst)+(tst*hst*hst/2))/Asp
                        Rsp = Rm-hst
                        Aeff = Asp*Rm/Rsp
                        pstar = 2*E*h*h/(Rm*Rm*np.sqrt(3*(1-v*v)))
                        V = P/pstar
                        n1 = 1 / 2 * np.sqrt(1 - V)
                        n2 = 1 / 2 * np.sqrt(1 + V)
                        L = Lsp - tst
                        Leff = 2*np.sqrt(Rm*h)/(3*(1-v*v))**(1/4)
                        degree = 2*L/Leff

                        F2 = ((np.cosh(n1 * degree) * np.sin(n2 * degree) / n2) + (

                        np.sinh(n1 * degree) * np.cos(n2 * degree) / n1)) / (

                        (np.cosh(n1 * degree) * np.sinh( n1 * degree) / n1) + (

                        np.cos(n2 * degree) * np.sin(n2 * degree) / n2))


                        F1 = 4/degree*(((np.cosh(n1*degree))**2-(np.cos(n2*degree))**2)/((np.cosh(n1*degree))*(np.sinh(n1*degree))/n1+(np.cos(n2*degree))*(np.sin(n2*degree))/n2))
                       
                        Wm = -(P*(Rm**2)/(E*h))*(1-(v/2))*(1-(Aeff*F2/(Aeff+tst*h+L*h*F1)))
                        ohmpim=E*(Wm/Rm)+(v*ohmx)

                        f = ohmx/ohmpim

                        ohmv = HY80

                        Et = E*(1-((ohmv-0.8*HY80orm)/(0.2*HY80orm))**2)
                        Es = E*ohmv/HY80orm/(0.8+0.2*math.atanh(np.tanh(0.9)))

                        vpbt = 0.5-(0.5-v)*Es/E

                        btz= ((2*(pi**2)*E*f)/(3*pii*(1-(v**2))))*((h/Rm)**2)*(((Rm*h)/(L**2))/(3-2*pii*(1-f)))    # 비대칭 탄성 좌굴강도

                        bbtz=btz*(1-v**2/1-vpbt**2)*(Et/E*(1-3*pii/4)+Es/E*3*pii/4)                                # 비대칭 비탄성 좌굴강도

                        k=ohmpim/ohmx

                        Kzegop  = 1-k+k**2

                        vp= 1/2-(Es/E)*(1/2-v)

                        H= 1+((1-(Et/Es))/(4*(1-vp**2)*Kzegop))*(((2-vp)-(1-2*vp)*k)**2-3*(1-vp**2))

                        A1 = 1 - ((1 - (Et  / Es)) /(4 * ( 1 - vp ** 2 ))) * Kzegop * H * (( 2 - vp ) - ( 1 - 2 * vp ) * k ) ** 2

                        A2 = 1 - (((1-(Et/Es))*(((1-2*vp)-((2-vp)*k))**2))/(4*(1-(vp**2))*Kzegop*H))

                        A12=1 + (1-(Et/Es))/(4*(1-(0.5-(0.5-v)*Es/E)**2)*Kzegop)*H

                        C=np.sqrt(((A1*A2)-(vpbt**2)*(A12**2))/(1-(vpbt**2)))

                        a= (3*((A2/A1)-vp**2*A12**2/A1**2)/(Rm**2*h**2))**(1/4)

                        dtz = (2 / np.sqrt(3 * (1 - v**2))) * E * (h**2 / Rm**2) * (((2 * L / (pi * Leff))**2) + ((1 / 4) * ((pi * Leff / (2 * L))**2))) #대칭 탄성 좌굴강도

                        dbtz = (2 / np.sqrt(3 * (1 - v**2))) * Es * ((h**2) / (Rm**2)) * C * (((a * L / (pi))**2 + (1 / 4) * (pi / (a * L))**2))  #대칭 비탄성 좌굴강도

                        r = 1-(0.25*(e**(-0.5*(btz/bbtz)-1)))  #감소계수

                        bbp =  r * bbtz # 비대칭 붕괴압력
                        dbp =  r * dbtz # 대칭 붕괴압력

                        중량=((h * Lsp) + Asp) * mildo * 10**-6

                        # 리스트에 추가
                        if bbp  >= 5.5 and dbp  >= 5.5:           #비대칭 붕괴압력, 대칭 붕괴압력이 6.0보다 크거나 5.5보다 작으면 탈락
                          if bbp  < 6.0 and dbp < 6.0:  
                             results.append([P, Rm, h, Lsp, hst, tst, bg, tg, bbp, dbp, 중량])

                       
# 표 출력
print("P\tRm\th\tLsp\thst\ttst\tbg\ttg\t비대칭붕괴압력\t\t대칭붕괴압력\t\t\t중량")
for row in results:
    print("\t".join(map(str, row)))

print("Results Length:", len(results))  # 케이스 개수
