import dspace


f = ['X1. = V1 + b21*X2 - b11*X1 - b12*X1',
     'X2. = V2 + b12*X1 + b32*X3 - b21*X2 - b22*X2 - b23*X2',
     'X3. = V3 + b23*X2 - b32*X3 - b33*X3']

eq = dspace.Equations(f)

ds = dspace.DesignSpace(eq, resolve_cycles=True)