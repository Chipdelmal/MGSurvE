import unittest
import MGSurvE.auxiliary as aux


###############################################################################
# Tests on distance
###############################################################################
class TestDistance(unittest.TestCase):
    def test_Vincenty(self):
        pointA = (6.250578, 53.478612)
        pointB = (5.916981, 50.752342)

        distanceM = aux.vincentyDistance(pointA, pointB, meters=True)
        testDistance = (304000 <= distanceM <= 304500)
        
        self.assertTrue(testDistance, msg='Expect Distance: {0}, Actual Distance: {1}'.format(304001, distanceM))

    def test_CheapRuler(self):
        pointA = (6.250578, 53.478612)
        pointB = (5.916981, 50.752342)

        distanceM = aux.cheapRuler(pointA, pointB)
        testDistance = (304000 <= distanceM <= 304500)
        
        self.assertTrue(testDistance, msg='Expect Distance: {0}, Actual Distance: {1}'.format(304001, distanceM))

    def test_Haversine(self):
        pointA = (6.250578, 53.478612)
        pointB = (5.916981, 50.752342)

        distanceM = aux.haversineDistance(pointA, pointB)
        testDistance = (304000 <= distanceM <= 304500)
        
        self.assertTrue(testDistance, msg='Expect Distance: {0}, Actual Distance: {1}'.format(304001, distanceM))


    def test_Vincenty_180Meridian(self):
        pointA = (177.221232, -17.947826)
        pointB = (-179.779055, -16.603513)

        distanceM = aux.vincentyDistance(pointA, pointB, meters=True)
        testDistance = (350500 <= distanceM <= 352000)
        
        self.assertTrue(testDistance, msg='Expect Distance: {0}, Actual Distance: {1}'.format(351826, distanceM))
    
    def test_CheapRuler_180Meridian(self):
        pointA = (177.221232, -17.947826)
        pointB = (-179.779055, -16.603513)

        distanceM = aux.cheapRuler(pointA, pointB)
        testDistance = (350500 <= distanceM <= 352000)
        
        self.assertTrue(testDistance, msg='Expect Distance: {0}, Actual Distance: {1}'.format(351826, distanceM))
    
    def test_Haversine_180Meridian(self):
        pointA = (177.221232, -17.947826)
        pointB = (-179.779055, -16.603513)

        distanceM = aux.haversineDistance(pointA, pointB)
        testDistance = (350500 <= distanceM <= 352000)
        
        self.assertTrue(testDistance, msg='Expect Distance: {0}, Actual Distance: {1}'.format(351826, distanceM))

    
    def test_Vincenty_0Distance(self):
        pointA = (177.221232, -17.947826)
        pointB = (177.221232, -17.947826)

        distanceM = aux.vincentyDistance(pointA, pointB, meters=True)
        testDistance = (distanceM == 0)
        
        self.assertTrue(testDistance, msg='Expect Distance: {0}, Actual Distance: {1}'.format(0, distanceM))

    def test_CheapRuler_0Distance(self):
        pointA = (177.221232, -17.947826)
        pointB = (177.221232, -17.947826)

        distanceM = aux.cheapRuler(pointA, pointB)
        testDistance = (distanceM == 0)
        
        self.assertTrue(testDistance, msg='Expect Distance: {0}, Actual Distance: {1}'.format(0, distanceM))

    def test_Haversine_0Distance(self):
        pointA = (177.221232, -17.947826)
        pointB = (177.221232, -17.947826)

        distanceM = aux.haversineDistance(pointA, pointB)
        testDistance = (distanceM == 0)
        
        self.assertTrue(testDistance, msg='Expect Distance: {0}, Actual Distance: {1}'.format(0, distanceM))
###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()