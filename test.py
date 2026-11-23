"""Small exact geometry regressions. No search output is needed."""
from copy import deepcopy
import unittest

from tile_cusps import decide
from verify import IDENTITY, act, cross, supporting_planes, verify_witness


class OctahedralTilings(unittest.TestCase):
    base=[(1,1,0,0),(1,-1,0,0),(1,0,1,0),(1,0,-1,0),(1,0,0,1),(1,0,0,-1)]

    def test_regular_octahedron_and_face_double(self):
        self.assertEqual(len(decide(self.base)["tiles"]),1)
        double=sorted(set(self.base)|{act(cross(IDENTITY,2),v) for v in self.base})
        witness=decide(double)
        self.assertEqual(len(witness["tiles"]),2)
        self.assertEqual(witness["six_times_klein_volume"],"248/27")
        # Opposite sides of the shared facet are necessary even when the two
        # original octahedra have the correct total volume.
        bad=deepcopy(witness)
        for malformed in ([1,0,0,0], [0,0,0,0,1]):
            bad["separators"]=[malformed]
            with self.assertRaises(ValueError):
                verify_witness(double,bad,supporting_planes(double))
        bad=deepcopy(witness)
        bad["tiles"][0]=[-1,*bad["tiles"][0][1:]]
        with self.assertRaises(ValueError):
            verify_witness(double,bad,supporting_planes(double))

    def test_octahedron_inside_a_larger_untileable_hull(self):
        result=decide(self.base+[(5,3,4,0)])
        self.assertEqual(result["possible_octahedra"],1)
        self.assertEqual(result["status"],"no octahedral tiling")
        self.assertEqual(result["six_times_klein_volume"],"44/5")

    def test_ideal_tetrahedron(self):
        result=decide([self.base[i] for i in [0,2,4,5]])
        self.assertEqual((result["possible_octahedra"],result["status"]),(0,"no octahedral tiling"))

    def test_projective_rescaling_and_invalid_input(self):
        rescaled=[tuple((i+2)*x for x in v) for i,v in enumerate(self.base)]
        self.assertEqual(decide(rescaled)["tiles"],decide(self.base)["tiles"])
        for bad in (self.base+[(1,0,0,0)],self.base+[self.base[0]],self.base[:3]):
            with self.assertRaises(ValueError):
                decide(bad)


if __name__=="__main__":
    unittest.main()
