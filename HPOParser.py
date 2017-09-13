
import pronto

def getTermsAndSyns(hpoFN):
    ont = pronto.Ontology(hpoFN)
    ret = []
    for t in ont:
        #print t
        #print '\t', t.id, t.name
        #ret.append((t.id, t.name))
        #for s in t.synonyms:
        #    ret.append((t.id, s.desc))
        fullList = [t.name]
        for s in t.synonyms:
            fullList.append(s.desc)
        
        ret.append((t.id, '; '.join(fullList)))
        
    return ret
    
if __name__ == '__main__':
    hpoFN = '/Users/matt/data/HPO_dl/hp.obo'
    terms = getTermsAndSyns(hpoFN)
    for t, d in terms:
        print t, d
    
    print str(len(terms))+' total elements'