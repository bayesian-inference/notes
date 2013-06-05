from numpy  import *

def get_new_idx(idx, card, old_stride, new_stride):
    new_idx = 0
    for i in range(len(old_stride)):
        new_idx += ((idx / old_stride[i]) % card[i]) * new_stride[i]
    #print "NEW_IDX=" + str(new_idx) + " OLD_STRIDE=" + str(old_stride) + " NEW_STRIDE=" + str(new_stride)
    return new_idx

class Factor:
    var = None
    card = None
    stride = None
    phi = None

    def __init__(self, var, card, phi):
        self.var = array(var)
        self.card = array(card)
        self.stride = concatenate([[1],cumprod(card)])
        self.phi = array(phi)

    def __str__(self):
        text = ""
        text += "VAR: " + str(self.var) + "\n"
        text += "CARD: " + str(self.card) + "\n"
        text += "VAL: " + str(self.phi)
        return text

    def __repr__(self):
        return str(self.var) + ", " + str(self.card) + ", " + str(self.phi)

    def margin(self, var):
        new_var = setdiff1d(self.var, var, assume_unique=True)
        old_idx = in1d(self.var,new_var).nonzero()

        #print "NEW_VAR " + str(new_var)
        #print "OLD_IDX " + str(old_idx)

        new_card = self.card[old_idx]
        new_stride = cumprod(concatenate([[1],new_card]))
        
        #print "NEW_CARD " + str(new_card)
        #print "NEW_STRIDE " + str(new_stride)

        new_phi = zeros(prod(new_card))
        for i in range(prod(self.card)):
            j = get_new_idx(i, new_card, self.stride[old_idx], new_stride)
            new_phi[j] += self.phi[i]

        return Factor(new_var, new_card, new_phi)

    def __mul__(self, fac2):
        if (self.var is None):
            return fac2
        if (fac2.var is None):
            return self
        var, index, inverse = unique(concatenate([self.var,fac2.var]), return_index=True, return_inverse=True)
        card = [self.card[i] if i < len(self.var) else fac2.card[i-len(self.var)] for i in index]
        inverse1 = inverse[0:len(self.var)]
        inverse2 = inverse[len(self.var):]

        psi = zeros(prod(card))

        j = 0
        k = 0
        assign = zeros(len(var))
        for i in range(prod(card)):
            #print "I-J-K: " + str((i,j,k))
            psi[i] = self.phi[j] * fac2.phi[k]
            #print "PSI_" + str(i) + "=" + str(psi[i])
            for l in range(len(var)):
                assign[l] += 1
                
                stride1_arr = where(inverse1 == l)[0]
                stride1 = 0
                if len(stride1_arr) > 0:
                    stride1 = self.stride[stride1_arr[0]]
                stride2_arr = where(inverse2 == l)[0]
                stride2 = 0
                if len(stride2_arr) > 0:
                    stride2 = fac2.stride[stride2_arr[0]]
                
                #print "STRIDE1: " + str(stride1) + ", STRIDE2: " + str(stride2)

                if assign[l] == card[l]:
                    assign[l] = 0
                    j -= (card[l]-1) * stride1
                    k -= (card[l]-1) * stride2
                    #print "==: l=" + str(l) + ", j=" + str(j) + ", k=" + str(k)
                else:
                    j += stride1
                    k += stride2
                    #print "!=: l=" + str(l) + ", j=" + str(j) + ", k=" + str(k)
                    break
        return Factor(var, card, psi)



def skuska():
    #f1 = Factor([1,2,4],[2,3,2],range(12))
    #f2 = Factor([1,3,4],[2,4,2],range(16))
    f1 = Factor([1,2],[2,2],[1,2,3,4])
    f2 = Factor([1,3],[2,2],[10,20,30,40])
    f = f1 * f2
    return f
        
