from pyliftover import LiftOver

lo = LiftOver('./AtoC.over.chain')

q = 0
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)
q = 23
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)

q = 24
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)
q = 43
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)

q = 55
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)

q = 56
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)
q = 66
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)
q = 77
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)

q = 78
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)
q = 97
x = lo.convert_coordinate('Achr1', q)
print(q,"\t",x)
