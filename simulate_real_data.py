import pysam
import os
proportion = [50, 55, 60, 65, 70, 75,80 , 85, 90, 95, 99]
# nums = [8796566, 7119889, 6906225, 6006724, 7474362]
nums = []
for i in range(1, 6):
    print(i)
    nums.append(int(pysam.view('-c', f'real_bam/in{i}.bam')))
print(nums)
print('FINISH NUMS')
for i in [2]: #range(1,11):
    bam_fname1=f'real_bam/in{i}.bam'
    for j in [3]:# range(i+1, 11):
        bam_fname2=f'real_bam/in{j}.bam'
        for v in [1]:
            for cov in [5, 10, 30]:
                n1 = nums[i-1]
                n2 = nums[j-1]
                for prop in proportion:
                    print(f'{cov},{i},{j},{prop},{v}')
                    print(n1,n2)
                    os.system(f'samtools view -h -b -s  {16569*cov*prop/(30*100*n1)} {bam_fname1} > 1.bam')
                    os.system(f'samtools view -h -b -s  {16569*cov*(1-prop)/(30*100*n2)} {bam_fname1} > 2.bam')
                    os.system(f'samtools cat 1.bam 2.bam | samtools sort > {i}-{j}_cov{cov}_{prop}_v{v}.bam')
            