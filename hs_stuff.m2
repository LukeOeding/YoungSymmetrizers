
-- hilbert series computations 
restart
R = ZZ[t] 
toDegRing = map(degreesRing R, R, {t=>T})
Q = (1-t^4)^5 * (1-t^6)*(1-t^8)^5*(1-t^10)*(1-t^12)^5
-- factor Q
-- Qp = (1-t^4)^5 *(1-t^6)*(1-t^8)^5*(1-t^12)^6
-- factor Qp
P=t^104+16*t^96+9*t^94+82*t^92 +145*t^90+383*t^88+770*t^86 +1659*t^84+3024*t^82 +5604*t^80+9664*t^78+15594*t^76 +24659*t^74 +36611*t^72+52409*t^70 +71847*t^68 +95014*t^66+119947*t^64+146849*t^62 +172742*t^60+195358*t^58+214238*t^56 +225699*t^54  +229752*t^52 +225699*t^50 +214238*t^48 +195358*t^46+172742*t^44 +146849*t^42+119947*t^40+95014*t^38 +71847*t^36+52409*t^34+36611*t^32 +24659*t^30 +15594*t^28 +9664*t^26+5604*t^24 +3024*t^22 +1659*t^20+770*t^18+383*t^16 +145*t^14+82*t^12+9*t^10+16*t^8+1    
reduceHilbert ( Divide(P,factor Q))
reduceHilbert ( Divide(toDegRing P,factor toDegRing Q))
factor P
peek dr 

peek Divide{P,Q}
factor(P/Q)
factor(P/Qp)
df = factor(Qp/Q)
factor(P/Qp)
P/Q - P/Q *Q/Qp * Qp/Q
factor(P*(Q/Qp))