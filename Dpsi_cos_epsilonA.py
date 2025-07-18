import math

def Dpsi_cos_epsilonA(T, F):
    """
    Calculate the approximate value of Delta psi * cos(epsilon_A) at TT Julian century
    T = (JD - 2451545)/36525

    F[] (input) are the fundamental arguments:
    F[0] = l = mean anomaly of the Moon
    F[1] = l' = mean anomaly of the Sun
    F[2] = L - Omega = mean longitude of the Moon - mean longitude of Moon's ascending node
    F[3] = D = mean elongation of the Moon from the Sun
    F[4] = Omega = mean longitude of Moon's ascending node
    F[5] = L_Me = mean longitide of Mercury
    F[6] = L_Ve = mean longitide of Venus
    F[7] = L_Ea = mean longitide of Earth
    F[8] = L_Ma = mean longitide of Mars
    F[9] = L_Ju = mean longitide of Jupiter
    F[10] = L_Sa = mean longitide of Saturn
    F[11] = L_Ua = mean longitide of Uranus
    F[12] = L_Ne = mean longitide of Neptune
    F[13] = p_A = general precession in longitude

    Delta psi is computed by a truncated IAU 2000A series, keeping terms with
    amplitudes greater than 6e-10 rad and 1e-11 rad/cy (122 terms).

    Returned Delta psi * cos(epsilon_A) in radians.
    """
    epsA = 0.4090926006005829 + T*(-0.00022707106390167 + T*(-8.876938501115605e-10 + T*(9.712757287348442e-09 + T*(-2.792526803190927e-12 - T*2.104091376015386e-13))))
    s = 8.341910002468069e-05*math.sin(F[4]+3.141398621411859)
    s += 6.38544187961502e-06*math.sin(2*F[2]-2*F[3]+2*F[4]-3.14055278673354)
    s += 1.103639471272347e-06*math.sin(2*F[2]+2*F[4]+3.140364408833701)
    s += 1.005772218329172e-06*math.sin(2*F[4]-0.0003364578476439898)
    s += 7.155482964182905e-07*math.sin(F[1]+0.008006593697976428)
    s += 3.44779871829835e-07*math.sin(F[0]-0.001226166810779262)
    s += 2.505620202702747e-07*math.sin(F[1]+2*F[2]-2*F[3]+2*F[4]-3.140578763247785)
    s += 1.877693986989801e-07*math.sin(2*F[2]+F[4]+3.140611507409127)
    s += 1.461544069750348e-07*math.sin(F[0]+2*F[2]+2*F[4]+3.138885869339794)
    s += 1.046368658184346e-07*math.sin(F[1]-2*F[2]+2*F[3]-2*F[4]+3.141078357595216)
    s += 7.611482188502514e-08*math.sin(F[0]-2*F[3]-3.140522576672403)
    s += 6.216626582066176e-08*math.sin(2*F[2]-2*F[3]+F[4]+0.001411558250493397)
    s += 5.985364333756062e-08*math.sin(F[0]-2*F[2]-2*F[4]+3.141438753852638)
    s += 3.072709235108941e-08*math.sin(2*F[3]-0.00236671010798622)
    s += 3.059659421492706e-08*math.sin(F[0]+F[4]+0.0004278244074266593)
    s += 2.891680223806012e-08*math.sin(F[0]-2*F[2]-2*F[3]-2*F[4]+0.00249810864368704)
    s += 2.810770733076956e-08*math.sin(F[0]-F[4]-0.003259958094348328)
    s += 2.50227666794999e-08*math.sin(F[0]+2*F[2]+F[4]+3.139093288480492)
    s += 2.313628013568719e-08*math.sin(2*F[0]-2*F[3]-0.0003771845091653969)
    s += 2.224955934316166e-08*math.sin(2*F[0]-2*F[2]-F[4]+3.140917169413882)
    s += 1.869748133639419e-08*math.sin(2*F[2]+2*F[3]+2*F[4]+3.13749580361748)
    s += 1.574723317611884e-08*math.sin(2*F[1]-2*F[2]+2*F[3]-2*F[4]+3.141592653589793)
    s += 1.505165953603724e-08*math.sin(2*F[0]+2*F[2]+2*F[4]+3.137373133440903)
    s += 1.417745186935223e-08*math.sin(2*F[0]-0.002530514723072703)
    s += 1.38622775924428e-08*math.sin(F[0]+2*F[2]-2*F[3]+2*F[4]-3.497359492156411e-05)
    s += 1.255041255256392e-08*math.sin(2*F[2]-0.002549536717151993)
    s += 1.056069829628443e-08*math.sin(2*F[2]-2*F[3]+3.140995857994301)
    s += 9.910077641445381e-09*math.sin(F[0]-2*F[2]-F[4]+3.141103440772307)
    s += 8.099783621226679e-09*math.sin(2*F[1]-0.0005985514338772515)
    s += 7.65715120853957e-09*math.sin(2*F[1]+2*F[2]-2*F[3]+2*F[4]-3.140579611008646)
    s += 7.351716594611941e-09*math.sin(F[0]-2*F[3]-F[4]+3.1408672514485)
    s += 6.81319431374123e-09*math.sin(F[1]+F[4]+3.13597113734224)
    s += 6.241032296017232e-09*math.sin(F[0]-2*F[3]+F[4]-3.138718428614093)
    s += 6.13490835259097e-09*math.sin(F[1]-F[4]+0.004978621738499853)
    s += 5.344590330395645e-09*math.sin(2*F[0]-2*F[2]-0.001269955775912256)
    s += 4.947053649557188e-09*math.sin(F[0]-2*F[2]-2*F[3]-F[4]+0.002450014698015139)
    s += 2.899014947660038e-09*math.sin(F[2]-F[3]+F[4]-F[7]-2*F[9]+5*F[10]+2.112678527143739)
    s += 3.726823817267823e-09*math.sin(F[0]+2*F[2]+2*F[3]+2*F[4]+3.135868766575786)
    s += 3.668104187987967e-09*math.sin(F[1]+2*F[2]+2*F[4]-0.001453871563522291)
    s += 3.563382666907906e-09*math.sin(F[0]+F[1]-2*F[3]-3.140504218645463)
    s += 3.462056669332664e-09*math.sin(F[1]-2*F[2]-2*F[4]+0.001120290807057238)
    s += 3.217731228665462e-09*math.sin(2*F[2]+2*F[3]+F[4]+3.137825909313233)
    s += 3.187671189169734e-09*math.sin(F[0]+2*F[3]-0.003650173902623447)
    s += 3.123656390929276e-09*math.sin(2*F[0]+2*F[2]-2*F[3]+2*F[4]-0.001086449983826601)
    s += 3.055295972212552e-09*math.sin(2*F[3]+F[4]+3.141275294032057)
    s += 2.811919517612434e-09*math.sin(F[0]+2*F[2]-2*F[3]+F[4]+0.0003448275725395338)
    s += 2.799323640768563e-09*math.sin(2*F[0]-2*F[3]-F[4]-0.002597846597843798)
    s += 2.593773175432264e-09*math.sin(2*F[0]+2*F[2]+F[4]+3.137667440104039)
    s += 2.395001224546236e-09*math.sin(2*F[3]-F[4]-0.004250986539198727)
    s += 2.303835071736334e-09*math.sin(F[1]-2*F[2]+2*F[3]-F[4]-0.0006313130474418824)
    s += 2.290746490151074e-09*math.sin(F[0]-F[1]-0.001269840587303581)
    s += 1.959348266025694e-09*math.sin(F[0]-F[3]-3.054136232211202)
    s += 2.107975460591217e-09*math.sin(F[1]-2*F[3]-3.13929274964128)
    s += 2.050763303757379e-09*math.sin(F[3]+3.140410621043379)
    s += 1.970769760485542e-09*math.sin(2*F[0]-2*F[3]+F[4]+0.001476013688257454)
    s += 1.966405784704713e-09*math.sin(F[0]-2*F[2]+0.001232740992911526)
    s += 1.361663502631758e-09*math.sin(2*F[9]-5*F[10]-F[13]+1.03080425520928)
    s += 1.735149857948379e-09*math.sin(F[1]+2*F[2]-2*F[3]+F[4]+0.00139703736997594)
    s += 1.643035353467799e-09*math.sin(F[0]+F[1]+3.140117293196702)
    s += 1.618805150347474e-09*math.sin(F[0]+2*F[2]-0.003893361579551144)
    s += 1.5882496933096e-09*math.sin(F[0]-F[1]-F[3]+3.141287403294024)
    s += 1.488863130424037e-09*math.sin(2*F[0]-2*F[2]-2*F[4]-0.0006512535712295179)
    s += 1.406463289717303e-09*math.sin(3*F[0]+2*F[2]+2*F[4]+3.136422068851762)
    s += 1.395299164784956e-09*math.sin(F[0]-F[1]+2*F[2]+2*F[4]+3.138812952618509)
    s += 1.3666939805692e-09*math.sin(F[0]+F[1]-2*F[2]-2*F[3]-2*F[4]+0.002483144949499867)
    s += 1.283312894781367e-09*math.sin(F[1]-2*F[2]-2*F[3]-2*F[4]+0.004155623981672819)
    s += 1.202827530382758e-09*math.sin(F[0]+F[1]+2*F[2]+2*F[4]-0.002821435479844587)
    s += 1.112173151408997e-09*math.sin(2*F[0]-F[4]-0.004359170295859138)
    s += 1.056409456124901e-09*math.sin(2*F[0]+F[4]-0.0009178519680427467)
    s += 8.092555506999434e-10*math.sin(4*F[7]-8*F[8]+3*F[9]+1.851235548038846)
    s += 1.001140251491192e-09*math.sin(3*F[6]-5*F[7]-2*F[13])
    s += 9.633291762249233e-10*math.sin(F[0]-2*F[2]+2*F[3]-F[4]-0.003019618401508975)
    s += 9.638095980457573e-10*math.sin(F[1]-F[2]+F[3]-F[4]-1.570796326794897)
    s += 9.604159022779907e-10*math.sin(F[0]+2*F[4]+3.141592653589793)
    s += 6.500127564502446e-10*math.sin(F[2]-F[3]+F[4]-8*F[6]+12*F[7]+0.4623133869773009)
    s += 8.047943613389794e-10*math.sin(2*F[2]+F[3]+2*F[4]-0.003012039083950873)
    s += 7.635870884552009e-10*math.sin(3*F[0]-0.003809505381148839)
    s += 7.374145180946386e-10*math.sin(F[0]-2*F[2]-4*F[3]-2*F[4]+0.005917090706027326)
    s += 7.199483164476608e-10*math.sin(F[6]-F[7])
    s += 6.981317007977317e-10*math.sin(8*F[7]-16*F[8]+4*F[9]+5*F[10])
    s += 6.811659824582501e-10*math.sin(F[0]-2*F[4]+3.138745686192619)
    s += 6.680739562168736e-10*math.sin(2*F[0]-2*F[2]-2*F[3]-2*F[4]-3.140141275799033)
    s += 6.48685234584542e-10*math.sin(F[0]-4*F[3]-3.137855750207207)
    s += 6.452986653777402e-10*math.sin(F[0]+2*F[2]+2*F[3]+F[4]+3.135582207560339)
    s += 6.370451769779302e-10*math.sin(F[0]-F[1]-F[3]-F[4]+3.141592653589793)
    s += 6.254096486313013e-10*math.sin(F[0]+F[1]+2*F[2]-2*F[3]+2*F[4])
    s += 6.215328409441865e-10*math.sin(2*F[0]-4*F[3]-3.139252564257516)
    s += 6.220159528635346e-10*math.sin(2*F[1]-2*F[2]+2*F[3]-F[4])
    s += 5.93093800438832e-10*math.sin(F[7]-F[9]-3.1178849131326)
    s += 6.050474740247008e-10*math.sin(2*F[2]-2*F[3]+3*F[4])
    s += 5.885688007494432e-10*math.sin(2*F[0]-2*F[2]-4*F[3]-2*F[4]+0.004118592857185193)
    s += T*8.444882361015366e-08*math.sin(F[4]+3.14142674108978)
    s += T*1.76333083538632e-09*math.sin(F[1]-3.137468512115215)
    s += T*7.943400401328656e-10*math.sin(2*F[2]-2*F[3]+2*F[4]+3.134268565802815)
    s += T*5.950611021781352e-10*math.sin(F[1]+2*F[2]-2*F[3]+2*F[4]+0.001629459206377831)
    s += T*2.397888466767765e-10*math.sin(F[1]-2*F[2]+2*F[3]-2*F[4])
    s += T*1.773939884115236e-10*math.sin(2*F[2]+F[4]+3.138859673239229)
    s += T*1.103963334611988e-10*math.sin(2*F[2]+2*F[4]+3.132809392415949)
    s += T*9.754451263923862e-11*math.sin(2*F[4])
    s += T*6.62255488395626e-11*math.sin(2*F[2]-2*F[3]+F[4])
    s += T*4.120916289431055e-11*math.sin(2*F[1]+3.141592653589793)
    s += T*3.490658503988658e-11*math.sin(2*F[1]+2*F[2]-2*F[3]+2*F[4])
    s += T*3.442177135877704e-11*math.sin(F[0])
    s += T*3.044629917367885e-11*math.sin(F[0]-F[4])
    s += T*3.044629917367885e-11*math.sin(F[0]+F[4])
    s += T*2.419220268736584e-11*math.sin(2*F[0]-2*F[2]-F[4]+3.141592653589793)
    s += T*2.031369323848956e-11*math.sin(F[0]+2*F[2]+F[4]+3.141592653589793)
    s += T*1.706544157505567e-11*math.sin(F[0]+2*F[2]+2*F[4]+3.141592653589793)
    s += T*1.21203420277384e-11*math.sin(F[1]+F[4]+3.141592653589793)
    s += T*1.018108730330026e-11*math.sin(F[1]+2*F[2]+2*F[4]+3.141592653589793)
    s += T*1.018108730330026e-11*math.sin(F[1]-2*F[2]-2*F[4]+3.141592653589793)
    s += T*1.01326059351893e-11*math.sin(F[0]-2*F[2]-F[4]+3.141592653589793)
    return s*math.cos(epsA)
