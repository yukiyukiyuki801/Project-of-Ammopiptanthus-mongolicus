#!/usr/bin/env python3
import importlib.util,tempfile,unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location('mapped_extent',Path(__file__).with_name('mapped_extent.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class ExtentTests(unittest.TestCase):
 def test_deletion_and_reference_skip_not_bridged(self):
  self.assertEqual(list(m.blocks(10,'5M3D4M20N2M')),[(10,15),(18,22),(42,44)])
 def test_insertions_and_clipping_do_not_inflate_reference(self):
  self.assertEqual(list(m.blocks(7,'3H5S4M2I3=1X2P4S')),[(7,15)])
 def test_invalid_cigar_rejected(self):
  for c in ('4MBOGUS','0M','*'):
   with self.assertRaises(AssertionError):list(m.blocks(0,c))
 def test_overlap_and_touch_merged_but_real_gap_retained(self):
  x=[('c',0,10),('c',4,15),('c',15,17),('c',18,20),('d',0,2)]
  self.assertEqual(list(m.merge_intervals(x)),[('c',0,17),('c',18,20),('d',0,2)])
 def test_unsorted_disjoint_rejected(self):
  with self.assertRaises(AssertionError):list(m.merge_intervals([('d',5,7),('c',0,2)]))
 def test_target_intersection_half_open_and_unique_length(self):
  rows=list(m.merge_intervals([('c',0,5),('c',3,10),('c',20,30)]));targets={'c':[3,5,10]}
  kept=list(m.select_regions(rows,targets));self.assertEqual(kept,[('c',0,10)]);self.assertEqual(sum(e-s for _,s,e in kept),10)
  covered={p for p in targets['c'] if any(c=='c' and s<=p<e for c,s,e in kept)}
  self.assertEqual(covered,{3,5});self.assertEqual(set(targets['c'])-covered,{10})
 def test_group_union_reuses_sample_beds_no_gap_fill(self):
  import heapq
  a=[('c',0,5),('c',20,25)];b=[('c',3,7),('c',26,28)]
  self.assertEqual(list(m.merge_intervals(heapq.merge(a,b))),[('c',0,7),('c',20,25),('c',26,28)])
 def test_fai_line_breaks_ambiguous_and_lowercase(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'ref.fa';p.write_bytes(b'>c\nACNT\ngtaC\n>other\nNN\n');Path(str(p)+'.fai').write_text('c\t8\t3\t4\t5\nother\t2\t20\t2\t3\n')
   f=m.Fasta(p)
   try:
    self.assertEqual(f.sequence('c',2,7),b'NTGTA');self.assertEqual(f.acgt('c',0,8),7);self.assertEqual(f.sequence('c',3,4),b'T');self.assertEqual(f.acgt('other',0,2),0)
    with self.assertRaises(AssertionError):f.sequence('c',7,9)
   finally:f.close()
if __name__=='__main__':unittest.main()
