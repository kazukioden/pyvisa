import cv2
import numpy as np
import matplotlib.pyplot as plt

def distance(pt1, pt2):
    return np.sqrt((pt1[0] - pt2[0])**2 + (pt1[1] - pt2[1])**2)
"""
start_frame = 11820
end_frame = 14790
path = r"E:\20221208_182653.avi"
"""
start_frame = 0
end_frame = 3000
path = r"C:\Users\anpan\OneDrive - The University of Tokyo\HD of Keyence\2023upto9monconference\20230607_154318.avi"   
cap = cv2.VideoCapture(path)
cap.set(cv2.CAP_PROP_POS_FRAMES, start_frame)
"""
points_roi1 = np.array([[325, 180], [345, 260], [415, 240], [395,160]])#roi1 = (320, 180, 70, 80)  # (x, y, width, height)roi2 = (390, 350, 70, 80)
points_roi2 = np.array([[370, 350], [390, 430], [460, 410], [440, 330]])
"""
points_roi1 = np.array([[800, 50], [800, 250], [1800, 250], [1800,50]])
points_roi2 = np.array([[650, 1100], [650, 1300], [1500, 1500], [1500, 1300]])

ft_params = dict(maxCorners=1000, qualityLevel=0.25, minDistance=20, blockSize=7)
lk_params = dict(winSize=(15, 15), maxLevel=2, criteria=(cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 1000, 0.001))

ret, frame = cap.read()
gray1 = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

def create_custom_rect_mask(image, points):
    mask = np.zeros_like(image, dtype=np.uint8)
    cv2.fillPoly(mask, [points], (255, 255, 255))
    return mask
mask_roi1 = create_custom_rect_mask(gray1, points_roi1)
mask_roi2 = create_custom_rect_mask(gray1, points_roi2)

ft1_roi1 = cv2.goodFeaturesToTrack(gray1, mask=mask_roi1, **ft_params)
ft1_roi2 = cv2.goodFeaturesToTrack(gray1, mask=mask_roi2, **ft_params)

print("ft1_roi1:", ft1_roi1)
print("ft1_roi2:", ft1_roi2)

mask = np.zeros_like(frame)



frame_count = 0
avg_distances = []

n = 0
while n < end_frame - start_frame:
    ret, frame = cap.read()
    if not ret:
        break
    frame_count += 1
    gray2 = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

    if ft1_roi1 is not None and ft1_roi2 is not None:
        ft2_roi1, status_roi1, err_roi1 = cv2.calcOpticalFlowPyrLK(gray1, gray2, ft1_roi1, None, **lk_params)
        ft2_roi2, status_roi2, err_roi2 = cv2.calcOpticalFlowPyrLK(gray1, gray2, ft1_roi2, None, **lk_params)
    else:
        break

    good1_roi1 = ft1_roi1[status_roi1 == 1]
    good2_roi1 = ft2_roi1[status_roi1 == 1]
    good1_roi2 = ft1_roi2[status_roi2 == 1]
    good2_roi2 = ft2_roi2[status_roi2 == 1]
    for i, (pt1_roi1, pt2_roi1) in enumerate(zip(good1_roi1, good2_roi1)):
        x1, y1 = pt1_roi1.astype(np.uint).ravel()
        x2, y2 = pt2_roi1.astype(np.uint).ravel()

        mask = cv2.line(mask, (x2, y2), (x1, y1), (0, 0, 200), 5)
        frame = cv2.circle(frame, (x2, y2), 5, (0, 0, 200), 3)

    for i, (pt1_roi2, pt2_roi2) in enumerate(zip(good1_roi2, good2_roi2)):
        x1, y1 = pt1_roi2.astype(np.uint).ravel()
        x2, y2 = pt2_roi2.astype(np.uint).ravel()

        mask = cv2.line(mask, (x2, y2), (x1, y1), (0, 0, 200), 5)
        frame = cv2.circle(frame, (x2, y2), 5, (0, 0, 200), 3)

    img = cv2.add(frame, mask)
   

    def overlay_custom_rectangle(image, points, color, alpha):
        overlay = image.copy()
        cv2.fillPoly(overlay, [points], color)
        return cv2.addWeighted(overlay, alpha, image, 1 - alpha, 0)

    color = (255, 255, 0)
    alpha = 0.3
    img = overlay_custom_rectangle(img, points_roi1, color=color, alpha=alpha)
    img = overlay_custom_rectangle(img, points_roi2, color=color, alpha=alpha)
    ratio = min(1000 / img.shape[1], 1000 / img.shape[0])
    new_size = (int(img.shape[1] * ratio), int(img.shape[0] * ratio))

    resized_img = cv2.resize(img, new_size)
    cv2.imshow('mask', resized_img)
    cv2.waitKey(1)
    if n==0:
        # 画像を保存する
        cv2.imwrite('new_image.png', resized_img)

    sorted_roi1 = sorted(good2_roi1, key=lambda pt: pt[0])
    sorted_roi2 = sorted(good2_roi2, key=lambda pt: pt[0])

    min_len = min(len(sorted_roi1), len(sorted_roi2))
    distances = [distance(sorted_roi1[i], sorted_roi2[i]) for i in range(min_len)]
    avg_distance = np.mean(distances)
    avg_distances.append(avg_distance)

    gray1 = gray2.copy()
    ft1_roi1 = good2_roi1.reshape(-1, 1, 2)
    ft1_roi2 = good2_roi2.reshape(-1, 1, 2)
    n += 1

    if cv2.waitKey(30) & 0xFF == ord('q'):#これ本当にいるか
        break

cv2.destroyAllWindows()
cap.release()

if len(avg_distances) > 0:
    time = np.arange(frame_count) / 15.0
    plt.plot(time, avg_distances/avg_distances[0])
    plt.xlabel("Time (s)")
    plt.ylabel("The Difference of Average Distance (%)")
    plt.show()
    path = path + 'OptFlow.txt'
    np.savetxt(path, np.c_[time, avg_distances], delimiter=',')
else:
    print("No average distances calculated.")

