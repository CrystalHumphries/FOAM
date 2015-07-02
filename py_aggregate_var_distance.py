import sys
import vcf

fh = sys.agv[1]

oldPos   = 0 
distance = 0
number   = 0

vcf_reader = vcf.Reader(filename=fh)
for record in vcf_reader.fetch('chr22', 0 , 90000000):
    if record.ALT != ".":
        number+=1
        if number > 1:
            dist = record.POS - oldPos
            distance+=dist
       oldPos = record.POS

print ("Sum:"      + str(distance))
avg_dist = distance/number
print ("Avg_dist:" + str(avg_dist))

    
        
