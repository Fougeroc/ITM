from sage.structure.sage_object import SageObject
from sage.plot.colors import rainbow
from sage.modules.free_module_element import free_module_element
from copy import copy
from sage.rings.integer import Integer
from sage.matrix.special import zero_matrix

def modulo(x,y,epsilon=0):
    if x >= 0:
        while x >= y-epsilon:
            x -= y
        return max(x,0)
    if x < 0:
        while x < -epsilon:
            x += y
        return max(x,0)

class IntervalTranslationMap(SageObject):
    def __init__(self, lengths, translations, labels=None, colors=None):
        if isinstance(lengths, dict) != isinstance(translations,dict):
            raise TypeError('lengths and translations should have the same format')
        if len(lengths) != len(translations) or (labels and len(labels) != len(lengths)) :
            raise TypeError('there should be the same number of lengths and translations')

        l, t = [], []
        if isinstance(lengths, dict):
            for lab in labels:
                l.append(lengths[lab])
                t.append(translations[lab])
        else:
            if not labels:
                labels = map(str,range(len(lengths)))

            for i in range(len(labels)):
                l.append(lengths[i])
                t.append(translations[i])

        try:
            self._ring = free_module_element(l+t).base_ring()
        except:
            raise TypeError('type of lengths and translation should be compatible')


        self._lengths = l 
        self._translations = t
        self._labels = labels
        if colors is None:
            self._colors = rainbow(len(self._labels), 'rgbtuple')
        else:
            self._colors = colors

        self.epsilon = 1e3*self._ring.epsilon()

    def __copy__(self):
        lengths = copy(self._lengths)
        translations = copy(self._translations)
        labels = copy(self._labels)
        colors = copy(self._colors)

        return(ITM(lengths, translations, labels, colors))

    def __len__(self):
        return len(self._labels)

    def vector(self, labels):
        v = []
        for l in labels:
            i = self._labels.index(l)
            v.append(self._lengths[i])
        for l in labels:
            i = self._labels.index(l)
            v.append(self._translations[i])
        return free_module_element(v)

    def plot_two_intervals(
        self,
        position=(0,0),
        vertical_alignment = 'center',
        horizontal_alignment = 'left',
        interval_height=1,
        labels_height=0.5,
        fontsize=14,
        labels=True,
        colors=None,
        latex=False,
        thickness=3,
        rauzy_move=False,
        alpha=.5,
        side=-1):
        r"""
        Returns a picture of the interval exchange transformation.

        INPUT:

        - ``position`` - a 2-uple of the position

        - ``horizontal_alignment`` - left (defaut), center or right

        - ``labels`` - boolean (defaut: True)

        - ``fontsize`` - the size of the label


        OUTPUT:

        2d plot -- a plot of the two intervals (domain and range)

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: t.plot_two_intervals()
            Graphics object consisting of 8 graphics primitives
        """
        from sage.plot.plot import Graphics, line2d, text

        G = Graphics()

        lengths = self._lengths
        total_length = sum(lengths)

        if colors is None:
            colors = self._colors

        s = self._ring.zero()
        if horizontal_alignment == 'left':
            s = position[0]
        elif horizontal_alignment == 'center':
            s = position[0] - total_length / 2
        elif horizontal_alignment == 'right':
            s = position[0] - total_length
        else:
            raise ValueError, "horizontal_alignement must be left, center or right"

        top_height = position[1] + interval_height
        for i in range(len(self._labels)):
            G += line2d([(s,top_height),(s+lengths[i],top_height)],
                        rgbcolor=colors[i],thickness=thickness, alpha=alpha)
            if labels == True:
                label_txt = str(self._labels[i]) 
                if latex:
                    label_txt = "$" + label_txt + "$"
                G += text(label_txt,
                          (s+lengths[i]/2,top_height+labels_height),
                          horizontal_alignment='center',
                          rgbcolor=colors[i],
                          fontsize=fontsize)

            s += lengths[i]

        if horizontal_alignment == 'left':
            s = position[0]
        elif horizontal_alignment == 'center':
            s = position[0] - total_length / 2
        elif horizontal_alignment == 'right':
            s = position[0] - total_length
        else:
            raise ValueError, "horizontal_alignement must be left, center or right"

        bottom_height = position[1] - interval_height
        ps = self._ring.zero()
        pos_left = s
        for i in range(len(self._labels)):
            s = pos_left + modulo(ps + self._translations[i], total_length, self.epsilon)
            if abs(self.total_length()-s) < self.epsilon: s=pos_left
            G += line2d([(s,bottom_height), (s+lengths[i],bottom_height)],
                        rgbcolor=colors[i],thickness=thickness, alpha=alpha)
            if labels == True:
                label_txt = str(self._labels[i])
                if latex:
                    label_txt = "$" + label_txt + "$"
                G += text(label_txt,
                          (s+lengths[i]/2,bottom_height-labels_height),
                          horizontal_alignment='center',
                          rgbcolor=colors[i],
                          fontsize=fontsize)
            ps += self._lengths[i]

        return G
    
    def plot_rauzy(self):
        position=[0,-5]
        plt = self.plot_two_intervals()
        I.rauzy(-1)
        plt += I.plot_two_intervals(position=(0,-5))

    def show(self, **kargs):
        r"""
        Shows a picture of the interval exchange transformation

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: phi = QQbar((sqrt(5)-1)/2)
            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,phi])
            sage: t.show()
        """
        normalized = copy(self)
        normalized.normalize()
        normalized.plot_two_intervals(**kargs).show(axes=False, figsize=(8,2))




    def total_length(self):
        return(sum(self._lengths))

    def length(self, lab):
        if isinstance(lab, int):
            return(self._lengths[lab])
        elif isinstance(lab, str):
            return(self._lengths[self._labels.index(lab)])
        else:
            raise TypeError('Label should be int or str')

    def lengths(self):
        return(self._lengths)

    def domain_singularities(self):
        s = [self._ring.zero()]
        for l in self._lengths:
            s.append(s[-1]+l)

    def in_which_interval(self, x):
        if (x < -self.epsilon or x >= self.total_length()+self.epsilon):
            raise ValueError("your value does not lie in [0; {}[".format(self.total_length()))
        if (abs(x - self.total_length()) < self.epsilon):
             return(False, len(self._labels))

        i = 0

        lengths = self._lengths
        
        while x >= -self.epsilon:
            x -= lengths[i]
            i += 1

        i -= 1
        x += lengths[i]

        if abs(x) <= self.epsilon:
            return(True, i)
        else:
            return(False, i)

    def rotation(self, c, inplace=True):
        from copy import copy
        singularity, i = self.in_which_interval(c)
        lengths = copy(self._lengths)
        translations = copy(self._translations)
        labels = copy(self._labels)

        if singularity:
            lengths = lengths[i:] + lengths[:i]
            translations = translations[i:] + translations[:i]
            labels = labels[i:] + labels[:i]
        else:
            l = c - sum(lengths[:i])
            t = translations[i]
            label = labels[i]

            lengths = [lengths[i]-l] + lengths[i+1:] + lengths[:i] + [l]
            translations = [t] + translations[i+1:] + translations[:i] + [t]
            labels = [label + 'r'] + labels[i+1:] + labels[:i] + [label + 'l']

        if inplace:
            self._labels = labels
            self._lengths = lengths
            self._translations = translations

        else:
            return ITM(lengths,translations,labels)
    
    def image(self, label):
        if isinstance(label, str): i = self._labels.index(label)
        elif isinstance(label, int): i = label
        else: raise TypeError('label ' + str(label) + ' should be an int or a str')
        l = self.total_length()
        s = sum(self._lengths[:i])
        left = modulo(s + self._translations[i], l)
        if abs(left) < self.epsilon: left = self._ring.zero()
        if abs(left-self.total_length()) < self.epsilon: left = self._ring.zero()
        right = modulo(s + self._translations[i] + self._lengths[i],l)
        if abs(right) < self.epsilon: right = self.total_length()
        if abs(right-self.total_length()) < self.epsilon: right = self.total_length()
        if right < 0: right += l
        return (left, right)

    def split(self, label):
        left, right = self.image(label)
        j = self._labels.index(label)

        if abs(self.total_length()-right) < self.epsilon\
        or abs(left) < self.epsilon:
            return None

        if right > left:
            raise ValueError('cannot split')

        l1, l2 = self.total_length()-left, right
        self._lengths[j] = l1
        self._lengths.insert(j+1, l2)
        self._labels.insert(j+1, self._labels[j] + 'r')
        self._colors.insert(j+1, self._colors[j])
        self._labels[j] += 'l'
        self._translations.insert(j, self._translations[j])

    def overlaping_intervals(self, side=-1):
        res = []
        for lab in self._labels:
            left, right = self.image(lab)
            if left > right or \
               side == -1 and abs(right-self.total_length()) < self.epsilon or \
               side == 0 and abs(left) < self.epsilon:
                   res.append(lab)
        return res

    def normalize(self):
        end = self.overlaping_intervals(1)
        for lab in end:
            self.split(lab)
        total_length = self.total_length()
        s = 0
        for i in range(len(self)):
            while total_length - self.epsilon < s + self._translations[i]: 
                self._translations[i] -= total_length
            s += self._lengths[i]
    
    def normalize_length(self, length=1):
        s = sum(self._lengths)
        x = length/s
        self._lengths = [l*x for l in self._lengths]
        self._translations = [t*x for t in self._translations]

    def symetric(self):
        self._labels = list(reversed(self._labels))
        self._lengths = list(reversed(self._lengths))
        self._colors = list(reversed(self._colors))
        self._translations = list(reversed(self._translations))
        self._translations = [-self._translations[i] for i in range(len(self))] 

    def is_rauzy(self, side=-1):
        labels = self._labels
        lengths = self._lengths
        translations = self._translations

        end = self.overlaping_intervals(side)
        if len(end) >= 2:
            return False
        if len(end) == 1:
            im = []
            for lab in self._labels:
                if not(end and  lab == end[0]):
                    left, right = self.image(lab)
                    im.append(right if side else left)

            next_s = max(im) if side else min(im)
            gap = self.total_length() - next_s if side else next_s
         
            lab = end[0]
            j = labels.index(lab)
            left, right = self.image(lab)
            l_top = lengths[side]
            l_bot = self.total_length() - left if side else right

            print min(l_bot,l_top), gap
            return(min(l_bot,l_top) < gap)

        if len(end) == 0:
            return True

    def rauzy(self, side=-1):
        self.normalize()
        end = self.overlaping_intervals(side)

        old_labels = copy(self._labels)
        labels = self._labels
        lengths = self._lengths
        translations = self._translations

        if len(end) >= 2:
            raise TypeError('Cannot apply rauzy induction in this case')

        im = []
        for lab in self._labels:
            if not(end and lab == end[0]):
                left, right = self.image(lab)
                im.append(right if side == -1 else left)
        next_s = max(im) if side == -1 else min(im)
        gap = self.total_length() - next_s if side else next_s

        if len(end) == 1:
            lab = end[0]
            j = labels.index(lab)
            left, right = self.image(lab)
            l_top = lengths[side]
            l_bot = (self.total_length() - left) if side == -1 else right

            if min(l_bot,l_top) > gap:
                raise ValueError('Overlaping does not unable induction')

            if abs(l_top-l_bot) < self.epsilon:
                translations[j] += self._translations[side]
                labels.pop(side)
                translations.pop(side)
                lengths.pop(side)

            if l_top > l_bot:
                lengths[side] -= l_bot
                translations[j] += translations[side]
                # labels[j] += '.' + self._labels[side]

            if l_top < l_bot:
                # create the preimage of the top interval
                lengths[j] -= l_top
                lengths.insert(j-side, l_top)
                labels.insert(j-side, labels[side])
                translations.insert(j-side, translations[j] + translations[side])
                # remove the top interval on the side
                labels.pop(side)
                lengths.pop(side)
                translations.pop(side)

        if len(end) == 0:
            if gap > self._lengths[side]:
                labels.pop(side)
                lengths.pop(side)
                translations.pop(side)
            else:
                lengths[side] -= gap

        colors = [self._colors[self._labels.index(lab)] for lab in old_labels]
        self._colors = colors

ITM = IntervalTranslationMap

def double_rotation(lengths,gap,p):
    from surface_dynamics.all import iet
    perm = iet.Permutation(p)
    l, t = dict(), dict()
    for i in range(len(perm)):
        l[perm[0][i]] = lengths[i]
        left = perm[1][0]
        i = 0
        t[left] = 0
        while perm[0][i] != left:
            t[left] -= lengths[i]
            i += 1
        right = perm[1][-1]
        i = 0
        t[right] = sum(lengths) - lengths[0]
        while perm[0][i] != right:
            i += 1
            t[right] -= lengths[i]
        middle = perm[1][1]
        t[middle] = gap
        i = 0
        while perm[0][i] != middle:
            t[middle] -= lengths[i]
            i += 1
    return ITM(l, t, perm[0])

def rauzy_type(itm):

    labels = itm._labels
    lengths = itm._lengths
    translations = itm._translations

    left, right = 0, 0
    for i in range(3):
        l, r = itm.image(i)
        if l < itm.epsilon: left = i
        if l > itm.total_length() - itm.epsilon: right = i

    middle = (set(range(3)) - {left,right}).pop()

    l_middle, r_middle = itm.image(middle)
    _, r_left = itm.image(left)

    if lengths[-1] < lengths[right]:
        return 'b0'
    if right == 0 and left == 2:
        if r_middle > r_left: 
            return 't1'
        elif l_middle > itm.length(left) - itm.length(middle) - itm.length(right):
            return 't2'
        else:
            return 't3'
    else:
        return 't3'

from surface_dynamics.interval_exchanges.labelled import LabelledPermutation
from surface_dynamics.all import iet

def interval_conversion(move):
    if move[0] == 'b':
        return 1, 0
    if move[0] == 't':
        x = int(move[1])
        return 0, x

class LabelledPermutationOver(SageObject):
    def __init__(self, p, labels=None):
        perm = iet.Permutation(p)
        self._perm = perm
        if labels: self._labels = labels
        else: self._labels = self._perm._alphabet

    def __repr__(self):
        return str(self._perm)

    def rauzy_matrix(self, move):
        from sage.matrix.special import identity_matrix
        n = len(self._perm)
        if move is None:
            return identity_matrix(2*n)

        side = -1

        winner, x = interval_conversion(move)

        ind = self._labels.rank
        winner_letter = ind(self._perm[winner][side])
        loser_letter = ind(self._perm[1-winner][side])
        free_letter = ind(self._perm[0][1])

        m = copy(identity_matrix(2*n))
        m[3+loser_letter, 3+winner_letter] = -1

        if x == 0:
            m[winner_letter, loser_letter] = 1
        if x == 1:
            m[winner_letter, winner_letter] = 0
            m[winner_letter, loser_letter] = 1
            m[winner_letter, loser_letter+3]  = 1
        if x == 2 or x == 3:
            m[winner_letter, winner_letter+3] = -1

        if move == 't2' or move == 't1':
            for i in range(n*2):
                for j in range(n, n*2):
                    m[i,j] = -m[i,j]

        return m

    def rauzy_move(self, move):
        perm = copy(self._perm)
        if move[0] == 'b':
            result = perm.rauzy_move('b')
        if move[0] == 't':
            if move[1] == '1':
                assert perm[0][0] == perm[1][-1] and perm[0][-1] == perm[1][0]
                result = (list(reversed(perm[0])), [perm[0][1], perm[0][0], perm[0][-1]])
            elif move[1] == '2':
                assert perm[0][0] == perm[1][-1] and perm[0][-1] == perm[1][0]
                result = (list(reversed(perm[0])), list(reversed(perm[1])))
            elif move[1] == '3':
                result = perm
        return LabelledPermutationOver(result, self._labels)
                        
    def moves(self):
        perm = self._perm
        if ((perm[0][0] == perm[1][-1]) and (perm[0][-1] == perm[1][0])):
                return ['b0','t1','t2','t3']
        else:
                return ['b0','t3']

    def random_move(self):
        import random
        return random.choice(self.moves())

    def right_matrix(self):
        ind = self._labels.rank
        perm = self._perm
        left_letter = perm[1][0]
        middle_letter = perm[1][1]
        right_letter = perm[1][2]
        m = zero_matrix(6,4)
        for i in range(3):
            m[i,i] = 1
        m[ind(middle_letter)+3,3] = 1
        for i in range(perm._twin[1][1]):
            m[ind(middle_letter)+3,ind(perm[0][i])] = -1
        for i in range(perm._twin[1][0]):
            m[ind(left_letter)+3,ind(perm[0][i])] = -1
        for i in range(perm._twin[1][2]+1,3):
            m[ind(right_letter)+3,ind(perm[0][i])] = 1
        return m

    def left_matrix(self):
        ind = self._labels.rank
        perm = self._perm
        m = zero_matrix(4,6)
        for i in range(3):
            m[i,i] = 1
        m[3,ind(self._perm[1][1])+3] = 1
        for i in range(perm._twin[1][1]):
            m[3,ind(perm[0][i])] = 1
        return m

    def alt_right_matrix(self):
        ind = self._labels.rank
        perm = self._perm
        left_letter = perm[1][0]
        middle_letter = perm[1][1]
        right_letter = perm[1][2]
        m = zero_matrix(6,4)
        for i in range(3):
            m[i,i] = 1
        m[2, 3] = 1
        m[ind(middle_letter)+3,3] = 1
        for i in range(perm._twin[1][1]):
            m[ind(middle_letter)+3,ind(perm[0][i])] = -1
            if perm[0][i] == left_letter:
                m[ind(middle_letter)+3,3] = -1
        for i in range(perm._twin[1][0]):
            m[ind(left_letter)+3,ind(perm[0][i])] = -1
            if perm[0][i] == left_letter:
                m[ind(left_letter)+3,3] = -1
        for i in range(perm._twin[1][2]+1,3):
            m[ind(right_letter)+3,ind(perm[0][i])] = 1
            if perm[0][i] == left_letter:
                m[ind(right_letter)+3,3] = 1
  
  
  
        return m

    def alt_left_matrix(self):
        ind = self._labels.rank
        perm = self._perm
        m = zero_matrix(4,6)
        for i in range(2):
            m[i,i] = 1
        m[2, ind(self._perm[1][0])] = 1
        m[3,ind(self._perm[1][1])+3] = 1
        m[2,ind(self._perm[1][1])+3] = -1
        for i in range(perm._twin[1][1]):
            m[3,ind(perm[0][i])] = 1
            m[2,ind(perm[0][i])] = -1
        return m



    def random(self):
        from sage.misc.prandom import random
        l = [random() for _ in range(3)]
        l.append(random()*l[self._perm[0].index(self._perm[1][0])])
        return l

    def matrix_path(self, s, restricted=True):
        from sage.matrix.special import identity_matrix
        R = self.left_matrix() if restricted else identity_matrix(6)
        aux = copy(self)
        n = len(s)
        for i in range(n/2):
            move = s[2*i:2*(i+1)]
            R = R*aux.rauzy_matrix(move)
            aux = aux.rauzy_move(move)
        if restricted:
            R = R*aux.right_matrix()
        return R

    def alt_matrix_path(self, s, restricted=True):
        from sage.matrix.special import identity_matrix
        R = self.alt_left_matrix() if restricted else identity_matrix(6)
        aux = copy(self)
        n = len(s)
        for i in range(n/2):
            move = s[2*i:2*(i+1)]
            R = R*aux.rauzy_matrix(move)
            aux = aux.rauzy_move(move)
        if restricted:
            R = R*aux.alt_right_matrix()
        return R


    def rauzy_path(self, s):
        aux = copy(self)
        n = len(s)
        for i in range(n/2):
            move = s[2*i:2*(i+1)]
            aux = aux.rauzy_move(move)
        return aux

    def __eq__(self, other):
        return(repr(self) == repr(other))

    def __hash__(self):
        return hash(repr(self))
