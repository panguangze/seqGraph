import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys

ref_list = {}
f = open(sys.argv[1])
for line in f:
    line = line.strip("\n").split()
    if line[1] not in ref_list:
        ref_list[line[1]] = int(line[4])
f.close()


ref_contig = {}
f_in = open(sys.argv[1])
title_contig = {}
ref_contig_l = {}


for line in f_in:
    # ref_length = ref_list[ref]
    contig_num = 10
    count = contig_num - 1
    contig = ""
    line = line.strip("\n").split("\t")
    if line[1] not in ref_contig.keys():
        ref_contig[line[1]] = []
        title_contig[line[1]] = []
        ref_contig_l[line[1]] = 0
    if line[0] == contig:
        start = min(int(line[10]), int(line[11]))
        stop = max(int(line[10]), int(line[11]))
        # rect = patches.Rectangle((start, count), width=(stop-start), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
        ref_contig[line[1]].append([start,stop,line[0]])
        if line[0] not in title_contig[line[1]]:
            title_contig[line[1]].append(line[0])
        ref_contig_l[line[1]] = ref_contig_l[line[1]] + (stop-start)
        # currentAxis.add_patch(rect)
        # title.add(line[0])
    else:
        count -= 1
        contig = line[0]
        start = min(int(line[10]), int(line[11]))
        stop = max(int(line[10]), int(line[11]))
        # rect = patches.Rectangle((start, count), width=(stop-start), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
        ref_contig[line[1]].append([start,stop,line[0]])
        if line[0] not in title_contig[line[1]]:
            title_contig[line[1]].append(line[0])
        ref_contig_l[line[1]] = ref_contig_l[line[1]] + (stop-start)
        # currentAxis.add_patch(rect)
        # title.add(line[0])
contig_ref = {}
# print(title_contig)

for ref in ref_list:
    title = set()
    ref_length = ref_list[ref]
    cover = [0]*ref_length
    vs = ref_contig[ref]
    # if ref == "CP024485.1":
    for v in vs:
        #print("CP024485.1",v[0],v[1])
        # if v[-1] in contig_cov.keys():
        #     contig_cov[v[-1]] = contig_cov[v[-1]] + (v[1]-v[0])
        # else:
        #     contig_cov[v[-1]] = (v[1]-v[0])
        for i in range(v[0],v[1]):
            cover[i] = 1
    # print(cover)
    un_covered = cover.count(0)
    # print(un_covered)
    # for c_k in contig_cov.keys():
    #     contig_l
    # contig_num = 10


    # f_in = open(sys.argv[1])
    # # plot
    # plt.figure(figsize=(20, 10))
    # plt.tight_layout()
    # plt.axis('off')
    # plt.xlim(xmin=0-5)
    # plt.xlim(xmax=ref_length+1000)
    # plt.ylim(ymin=-0.5)
    # plt.ylim(ymax=contig_num+2)

    # currentAxis = plt.gca()
    # rect = patches.Rectangle((1, contig_num), width=ref_length, height=0.3, linewidth=2, edgecolor="black", facecolor="#549FCF")
    # currentAxis.add_patch(rect)
    # contig = ""
    # count = contig_num - 1
    # print(un_covered)

    # un_covered = 0
    if ref in ref_contig.keys():
        # if (float(ref_contig_l[ref])/float(ref_length))<0.8:
        if un_covered / ref_length>0.4:
            print("not ok",ref,un_covered,un_covered / ref_length,ref_length)
            plt.close()
            f_in.close()
            continue



    #     vs = ref_contig[ref]
    #     print(ref,ref_contig_l[ref],ref_length,float(ref_contig_l[ref])/float(ref_length))
    #     for v in vs:
    #         if v[-1] == contig:
    #             rect = patches.Rectangle((v[0], count), width=(v[1]-v[0]), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
    #             currentAxis.add_patch(rect)
    #         else:
    #             count -= 1
    #             contig = v[-1]
    #             rect = patches.Rectangle((v[0], count), width=(v[1]-v[0]), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
    #             currentAxis.add_patch(rect)

    # for k in ref_contig.keys():
    #     vs = ref_contig[k]
    #     for v in vs:
    #         currentAxis.add_patch(rect)
    #     line = line.strip("\n").split("\t")
    #     if line[1] == ref:
    #         if line[0] == contig:
    #             start = min(int(line[10]), int(line[11]))
    #             stop = max(int(line[10]), int(line[11]))
    #             rect = patches.Rectangle((start, count), width=(stop-start), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
    #             currentAxis.add_patch(rect)
    #             title.add(line[0])
    #         else:
    #             count -= 1
    #             contig = line[0]
    #             start = min(int(line[10]), int(line[11]))
    #             stop = max(int(line[10]), int(line[11]))
    #             rect = patches.Rectangle((start, count), width=(stop-start), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
    #             currentAxis.add_patch(rect)
    #             title.add(line[0])
    pt=""
    for i in title_contig[ref]:
        pt = pt+i
    if pt in contig_ref.keys():
        contig_ref[pt].append(ref)
    else:
        contig_ref[pt] = [ref]
    # plt.title(pt)
    # plt.savefig("%s_blast.png"%ref, dpi=600)
    # plt.close()
    # f_in.close()
# for k in contig_ref.keys():
#     print(k,[ref_list[i] for i in contig_ref[k]], contig_ref)
import re
k_lens = {}
for k in contig_ref.keys():
    k_lens[k] = []
    t_l = 0
    t = re.split(r"[+-]",k.strip())
    # print(t,"ewwe")
    for i in t:
        if i == "":
            continue
        l = int(i.split("_")[3])
        k_lens[k].append(l)

# print("k_lens===============================================")
# print(k_lens)
# print("k_lens===============================================")

import math
count = 0
result = []
skip = []
replace = {}
for fk in k_lens.keys():
    if fk in skip:
        continue
    a = k_lens[fk]
    flag = True
    for sk in k_lens.keys():
        b = k_lens[sk]
        if fk == sk or sk < fk or sk in skip:
            continue
        tmp = [j for j in a if j in b]
        if sum(tmp)/sum(a) > 0.9 and sum(tmp)/sum(b) > 0.9:
            # min_diff_ref_with_contig = 1
            # min_ref = ""
            # for ref in contig_ref[fk]:
            #     if math.abs(1-ref_list[ref]/sum(a)) >
            # pass
            print(sum(tmp)/sum(a),sum(tmp)/sum(b))
            count = count+1
            # flag = False
            # print(fk,"----------",sk,"duplicate",tmp)
            # replace[fk] = sk
            replace[sk] = fk
            # if sum(b) < sum(a):
            #     if sk not in result:
            #         result.append(sk)
            #     skip.append(fk)
            # else:
            #     if fk not in result:
            #         result.append(fk)
            #     skip.append(sk)
        # else:

    if flag:
        result.append(fk)

# print(count)
# # print(tmp)
# print(len(contig_ref))
# print(result)

t = 0
for k in result:
    if "self" in k or "gene" in k:
        pass
        # print(k)
    t = t+ len(contig_ref[k])
    for ref in contig_ref[k]:
        # pass

        ref_length = ref_list[ref]
        contig_num = 30
        cover = list(range(0,ref_length))
        # plot
        plt.figure(figsize=(20, 10))
        plt.tight_layout()
        plt.axis('off')
        plt.xlim(xmin=0-5)
        plt.xlim(xmax=ref_length+1000)
        plt.ylim(ymin=-0.5)
        plt.ylim(ymax=contig_num+2)

        currentAxis = plt.gca()
        rect = patches.Rectangle((1, contig_num), width=ref_length, height=0.3, linewidth=2, edgecolor="black", facecolor="#549FCF")
        currentAxis.add_patch(rect)
        contig = ""
        count = contig_num - 1
        vs = ref_contig[ref]
        # if ref == "CP077390.1":
        # for v in vs:
        #     print(v[0],v[1])
        #     for i in range(v[0],v[1]):
        #         cover[i] = 1
        # un_covered = cover.count(0)
        # print(ref,ref_contig_l[ref],ref_length,float(ref_contig_l[ref])/float(ref_length))
        for v in vs:
            if (v[1]-v[0])/ref_length < 0.05 and (v[1]-v[0]) < 300:
                continue
            if v[-1] == contig:
                rect = patches.Rectangle((v[0], count), width=(v[1]-v[0]), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
                currentAxis.add_patch(rect)
            else:
                count -= 1
                contig = v[-1]
                rect = patches.Rectangle((v[0], count), width=(v[1]-v[0]), height=0.3, linewidth=2, edgecolor="black", facecolor="#FF8C00")
                currentAxis.add_patch(rect)
        #plt.title(pt)
        if k in replace.keys():
            k2 = replace[k]
        else:
            k2 = k
        print("res: ", k2)
        plt.savefig("N"+k2[0:15]+"_"+k2[-15:-1]+"%s_blast.png"%ref, dpi=600)
        plt.close()
        f_in.close()
print(t)