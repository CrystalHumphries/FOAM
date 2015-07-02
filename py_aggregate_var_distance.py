import sys
import vcf
import tabix


fh = sys.argv[1]

oldPos   = 0 
distance = 0
number   = 0

#vcf_reader = vcf.Reader(filename=fh)
tb  = tabix.open(fh)
records = tb.query('chr22', 0, 90000000)

for record in records:
    if record[4] is not '.':
        number+=1
        pos = int(record[1])
        if number > 1:
            dist = pos - oldPos
            distance+=dist
        oldPos = pos

print ("Sum:"      + str(distance))
avg_dist = distance/number
print ("Avg_dist:" + str(avg_dist))

    
        
