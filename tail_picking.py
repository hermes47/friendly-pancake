source = '/Users/iwelsh/Downloads/LMSDSearchResultsDownload13H10M37S07Dec15.csv'
with open(source,'r') as fh:
    data = fh.readlines()[1:]
tails = set()
lipids = []
for line in data:
    s_l = line.split('","')
    if not (s_l[1].startswith('PC(') or s_l[1].startswith('PE(') or s_l[1].startswith('PG(') or s_l[1].startswith('PS(') or s_l[1].startswith('PI(') or s_l[1].startswith('PA(')):
        continue
    lipids.append(s_l[1])
    com_name = s_l[1][3:-1]
    tails.update(com_name.split('/'))
tails = list(tails)
for i in sorted(tails):
    print(i)

with open('/Users/iwelsh/GitHub/itchy-wookie/Lipid_Input/needed_tails.txt','r') as fh:
    load_tails = fh.readlines()
count = 0
used = set()
for hg in ['PC','PE','PG','PS','PI','PA']:
    for t1 in load_tails:
        for t2 in load_tails:
            if '{}({}/{})'.format(hg,t1.strip(),t2.strip()) in lipids:
                count += 1
                used.update([t1.strip(),t2.strip()])
print(count,'out of ',len(load_tails)**2*6)
print('unused tails: ',set([x.strip() for x in load_tails]).difference(used))
