#!/usr/bin/env python

import re
def parse_pileup_bases(bases):
  base_pattern = re.compile(r"^(\^[^^])?[ATGCatgc,.*Nn]\$?")
  indel_pattern = re.compile(r"^(\+|\-)")
  indel_with_number_pattern = re.compile(r"^(\+|\-)[1-9][0-9]*")
  parsed = []
  s = bases
  last_base = ""
  while s != "":
    if indel_pattern.match(s):
      m = indel_with_number_pattern.search(s)
      len_indel = int(s[:m.end()][1:])
      parsed[-1] = last_base + s[:m.end() + len_indel]
      s = s[m.end() + len_indel:]
    else:
      m = base_pattern.search(s)
      last_base = s[:m.end()]
      parsed.append(last_base)
      s = s[m.end():]
  return parsed

def read_pileup(file):
  columns = [("sequence_identifier", str),
             ("position", int),
             ("reference_nucleotide", str),
             ("num_reads", int),
             ("bases", parse_pileup_bases),
             ("quolities", str)]
  for line in file:
    d = {}
    for col,data in zip(columns,line.split()):
      d[col[0]] = col[1](data)
    yield d

def count_pileup_base(data):
  read_begin_end_pattern = re.compile(r"(\^.)|\$")
  ref_match_pattern = re.compile(r"[.,]")
  ref = data["reference_nucleotide"]
  res = {}
  for b in data["bases"]:
    b = ref_match_pattern.sub(ref, read_begin_end_pattern.sub("",b.upper()))
    res[b] = res.get(b,0) + 1
  return res

def major_snps_indels(pileup_file):
  insertion = re.compile(r"^[ATGCatgc.,]\+[1-9][0-9]*[ATGCatgc]+")
  deletion = re.compile(r"^[ATGCatgc.,]-[1-9][0-9]*[ATGCatgc]+")
  d = {}
  for line in pileup_file:
    ref = line["reference_nucleotide"]
    position = line["position"]
    num_reads = line["num_reads"]
    count = count_pileup_base(line)
    major = sorted(count.items(),
                   key = lambda l: l[1],reverse=True)[0]
    major_nuc = major[0]
    num_major = major[1]
    if len(major_nuc) == 1 and major_nuc != ref:#snp
      d[position] = {"nucleotide": major_nuc,
                     "num_all_reads": num_reads,
                     "num_major_reads": num_major}
    elif insertion.match(major_nuc):#insertion
      nucleotide = re.sub(r"\+[1-9][0-9]*", "", major_nuc)
      d[position] = {"nucleotide": nucleotide,
                     "num_all_reads": num_reads,
                     "num_major_reads": num_major}
    elif deletion.match(major_nuc):#deletion
      deletion_length = int(re.findall(r"\d+",major_nuc)[0])
      nucleotide = major_nuc.split("-")[0].upper()
      if nucleotide != ref:
        d[position] = {"nucleotide": nucleotide,
                       "num_all_reads": num_reads,
                       "num_major_reads": num_major}
  return d

def read_sam(file):
  columns = [("qname", str),
             ("flag", int),
             ("rname", str),
             ("pos", int),
             ("mapq", int),
             ("cigar", parse_cigar),
             ("rnext", str),
             ("pnext", int),
             ("tlen", int),
             ("seq", str),
             ("qual", str)]
  for line in file:
    line = line.strip()
    if line[0] != "@":
      d = {}
      for col,data in zip(columns,line.split("\t")[:len(columns)]):
        d[col[0]] = col[1](data)
      d["option"] = "\t".join(line.split("\t")[len(columns):])
      yield d

def parse_cigar(cigar):
  pattern = re.compile(r"[1-9][0-9]*[MIDNSHP=X]")
  return map(lambda s: (int(s[:-1]), s[-1]), pattern.findall(cigar))

def arrange_seq(sam_read):
  seq = sam_read["seq"]
  pos = sam_read["pos"]
  cigar = sam_read["cigar"]
  clipped_seq = seq
  if cigar[0][1] == "S":
    len = cigar[0][0]
    clipped_seq = clipped_seq[len:]
  if cigar[-1][1] == "S":
    len = cigar[-1][0]
    clipped_seq = clipped_seq[:-len]
  res = []
  for c in cigar:
    len = c[0]
    if c[1] == "M":
      res += list(clipped_seq[:len])
      clipped_seq = clipped_seq[len:]
    elif c[1] == "I":
      res[-1] = res[-1] + clipped_seq[:len]
      clipped_seq = clipped_seq[len:]
    elif c[1] == "D":
      res += list(["*"] * len)
  return res

def filter_minor_reads_on_variations(sam,
                                     snp_indel,
                                     threshold,
                                     spec_positions=[],
                                     invert=False):
  d = {}
  for p in snp_indel.keys():
    num_major = snp_indel[p]["num_major_reads"]
    num_all = snp_indel[p]["num_all_reads"]
    if num_major/float(num_all) > threshold:
      d[p] = snp_indel[p]
  snp_indel = d

  for read in sam:
    if not read["flag"] & 4:
      begin_pos = read["pos"]
      seq = arrange_seq(read)
      end_pos = begin_pos + len(seq)
      included_snp_positions = [ p for p in snp_indel.keys() if p >= begin_pos and p < end_pos]
      if spec_positions:
        included_snp_positions = sorted(list(set(included_snp_positions) & set(spec_positions)))
      if invert:
        inv = lambda f: not f
      else:
        inv = lambda f:f
      f = True
      for p in included_snp_positions:
        sam_nuc = seq[p - begin_pos]
        if sam_nuc != snp_indel[p]["nucleotide"]:
          f = False
      if inv(f): 
        yield read

import sys
if __name__ == "__main___":
  for line in read_pileup(open(sys.argv[1])):
    bases = line["bases"]
    ref = line["reference_nucleotide"]
    print ref
    print line["num_reads"]
    print bases
    print count_pileup_base(line)
    print 
    num_check = line["num_reads"] != len(bases)
    if num_check:
      print "wrong", line, parsed

if __name__ == "__main___":
  m = major_snps_indels(read_pileup(open(sys.argv[1])))
  for p in sorted(m.keys()):
    print p, m[p]["nucleotide"], m[p]["num_all_reads"], m[p]["num_major_reads"]

if __name__ == "__main___":
  for line in read_sam(open(sys.argv[1])):
    print line["cigar"]
    print ",".join(arrange_seq(line))

import argparse
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("sam")
  parser.add_argument("pileup")
  parser.add_argument("-t", "--remove_threshold", default = 0.7, type=float)
  parser.add_argument("-v", "--invert", action="store_true", default=False)
  parser.add_argument("-p", "--position_file", type=argparse.FileType("r"))
  args = parser.parse_args()
  
  spec_pos = []
  if args.position_file:
    spec_pos = [ int(l.strip().split(":")[0]) for l in args.position_file if l ]
  
  for line in open(args.sam):
    if line[0] == "@":
      print line.strip()
    else:
      break
  
  sam_format = "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\t{option}"
  sam = read_sam(open(args.sam))
  snp_indel = major_snps_indels(read_pileup(open(args.pileup))) 
  for read in filter_minor_reads_on_variations(sam,
                                               snp_indel,
                                               args.remove_threshold,
                                               spec_pos,
                                               args.invert):
    read["cigar"] = "".join(str(t[0]) + t[1] for t in read["cigar"])
    print sam_format.format(**read)
