import cv2
import numpy as np

# 画像読み込み
img = cv2.imread("20230313_191353renzokukiridashi.jpg")#"20230306_163322.jpg""E:\20230306_171207.jpg""E:\20230306_173031.jpg""C:\Users\anpan\OneDrive - The University of Tokyo\デスクトップ\python_lesson\20230313_191316renzokukiridashiarea.jpg"

# RGB色空間に変換
rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

# ガウシアンフィルタを適用して平滑化
blur = cv2.GaussianBlur(rgb, (21, 21), 5)

# 薄いピンク色の範囲を指定
lower_pink = np.array([230, 0, 0])  #初期lower_pink = np.array([100, 0, 0])  upper_pink = np.array([255, 200, 255])

upper_pink = np.array([255, 180, 255])

# 指定した範囲内の色だけを抽出するマスク画像の生成
mask = cv2.inRange(blur, lower_pink, upper_pink)
cv2.imshow('result', mask)
cv2.waitKey(1000)
# 輪郭検出
contours, hierarchy = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

# 最も長い輪郭を取得
max_contour = max(contours, key=lambda cnt: cv2.arcLength(cnt, True))

# 輪郭の周長を計算
perimeter = cv2.arcLength(max_contour, True)

# 1000角形での近似
epsilon = 0.0001 * perimeter
approx = cv2.approxPolyDP(max_contour, epsilon, True)

# 面積を計算
area = cv2.contourArea(approx)

# 結果の表示
result = img.copy()
cv2.drawContours(result, [approx], 0, (0, 255, 0), 2)
cv2.putText(result, f"Area: {area:.2f} pixels", (500, 500), cv2.FONT_HERSHEY_SIMPLEX, 2, (0, 0, 0), 2)
cv2.imshow('result', result)
cv2.waitKey(0)
cv2.imwrite('20230313_191353renzokukiridashiarea.jpg', result)
cv2.destroyAllWindows()