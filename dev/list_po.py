import re, sys
s = open(sys.argv[1], encoding='utf-8').read()
ids = re.findall(r'msgid "((?:[^"\\]|\\.)*)"\nmsgstr', s)
ids = [i for i in ids if i]
print(len(ids))
for i in ids:
    print(repr(i))
