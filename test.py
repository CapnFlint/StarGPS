import numpy as np
from mlat import MLAT


def distance(a, b):
    return np.linalg.norm(a - b)

def start():
    markers = np.array([[0,0,0],[5,0,0]])
    # Marker 1, bottom left of triangle
    om1 = np.array([0, 0, 0])
    # Marker 2, bottom right of triangle
    om2 = np.array([5,0,0])
    # Start point, in the middle
    test = np.array([2.5,2.5,2.5])

    # Ranges to make a right angle triangle. Should give coordinates: 5, 0, 5
    ranges = np.array([5,7.07107])

    for i in range(10): 
        deltas = np.empty(3)
        for j in range(markers.shape[0]):
            deltas[j] = distance(markers[j], test)
            print test
            print markers[j]
            print distance(test, markers[j])
            print test - markers[j]
            print "-----"
        break

        for j in range(markers.shape[0]):
            delta = np.zeros(3)
            dist = distance(test, markers[j])
            print "To OM{0}: {1}".format(j, str(dist))
            adjust = (ranges[j] - dist) / dist * (test - markers[j])
            print "Adjustment"
            print adjust
            delta += adjust
            print delta
            delta *= 2 * 0.01
            print "Delta:"
            print delta
            test1 = test + delta
            if distance(markers[j], test1) < dist:
                test = test1
            else:
                test = test - delta


if __name__ == "__main__":
    start()


