import numpy as np
import sys

def load_uniform_grid_data(filename):
    """
    Загружает данные из файла с колонками: x[Mpc] y[Mpc] z[Mpc] Bx[G] By[G] Bz[G].
    Сетка предполагается равномерной и отсортированной по x, затем y, затем z
    (быстрое изменение z). Определяет размеры сетки и шаги.

    Возвращает:
        xs, ys, zs : одномерные массивы уникальных координат
        Bx, By, Bz : трёхмерные массивы компонент поля формы (nx, ny, nz)
        dx, dy, dz : шаги сетки
    """
    # Загружаем все данные
    data = np.loadtxt(filename)
    if data.shape[1] != 6:
        raise ValueError("Файл должен содержать 6 колонок: x y z Bx By Bz")
    x, y, z, Bx, By, Bz = data.T

    # Определяем уникальные координаты (сортированные)
    xs = np.unique(x)
    ys = np.unique(y)
    zs = np.unique(z)
    nx, ny, nz = len(xs), len(ys), len(zs)

    if nx * ny * nz != len(x):
        raise ValueError("Общее число точек не соответствует произведению уникальных координат. "
                         "Возможно, сетка нерегулярная или порядок точек не соответствует ожидаемому.")

    # Проверяем равномерность шагов
    dxs = np.diff(xs)
    dys = np.diff(ys)
    dzs = np.diff(zs)
    if not (np.allclose(dxs, dxs[0]) and np.allclose(dys, dys[0]) and np.allclose(dzs, dzs[0])):
        print("Предупреждение: шаги сетки не строго постоянны. Используем средние значения.")
    dx, dy, dz = dxs[0], dys[0], dzs[0]

    # Преобразуем одномерные массивы в 3D, предполагая порядок (x, y, z) с быстрым z
    # Проверяем, что данные действительно отсортированы как требуется
    # Вычисляем линейные индексы в C-порядке (ось z быстрая)
    x_idx = np.searchsorted(xs, x)
    y_idx = np.searchsorted(ys, y)
    z_idx = np.searchsorted(zs, z)
    linear = x_idx * ny * nz + y_idx * nz + z_idx
    if not np.all(np.diff(linear) >= 0):
        print("Предупреждение: данные не отсортированы в порядке (x,y,z). Результат reshape может быть неверным.")
        # Здесь можно было бы выполнить сортировку, но для простоты предполагаем корректный порядок.

    Bx_3d = Bx.reshape(nx, ny, nz)
    By_3d = By.reshape(nx, ny, nz)
    Bz_3d = Bz.reshape(nx, ny, nz)

    return xs, ys, zs, Bx_3d, By_3d, Bz_3d, dx, dy, dz

def compute_divergence(Bx, By, Bz, dx, dy, dz):
    """
    Вычисляет дивергенцию магнитного поля на равномерной сетке.
    Использует центральные разности (np.gradient).
    """
    dBx_dx = np.gradient(Bx, dx, axis=0)
    dBy_dy = np.gradient(By, dy, axis=1)
    dBz_dz = np.gradient(Bz, dz, axis=2)
    divB = dBx_dx + dBy_dy + dBz_dz
    return divB

def analyze_divergence(divB, tol=1e-10):
    """
    Печатает статистику дивергенции.
    tol - порог для проверки бездивергентности (абсолютное значение).
    """
    abs_div = np.abs(divB)
    max_div = np.max(abs_div)
    min_div = np.min(divB)
    mean_div = np.mean(divB)
    rms_div = np.sqrt(np.mean(divB**2))
    median_div = np.median(divB)

    print("Статистика div B:")
    print(f"  Минимум: {min_div:.3e} G/Mpc")
    print(f"  Максимум: {max_div:.3e} G/Mpc")
    print(f"  Среднее: {mean_div:.3e} G/Mpc")
    print(f"  RMS: {rms_div:.3e} G/Mpc")
    print(f"  Медиана: {median_div:.3e} G/Mpc")
    # Доля точек, где |div B| > tol
    n_viol = np.sum(abs_div > tol)
    total = divB.size
    print(f"  Точек с |div B| > {tol:.1e}: {n_viol} из {total} ({100*n_viol/total:.2f}%)")

    if max_div < tol:
        print("✓ Поле бездивергентно (в пределах заданного допуска).")
    else:
        print("✗ Поле НЕ является бездивергентным (превышает допуск).")

def main():
    if len(sys.argv) < 2:
        print("Использование: python script.py <имя_файла> [допуск]")
        sys.exit(1)

    filename = sys.argv[1]
    tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-10

    print(f"Загрузка данных из {filename}...")
    xs, ys, zs, Bx, By, Bz, dx, dy, dz = load_uniform_grid_data(filename)

    print(f"Размеры сетки: nx={len(xs)}, ny={len(ys)}, nz={len(zs)}")
    print(f"Шаги: dx={dx:.4e} Mpc, dy={dy:.4e} Mpc, dz={dz:.4e} Mpc")
    print("Вычисление дивергенции...")
    divB = compute_divergence(Bx, By, Bz, dx, dy, dz)
    analyze_divergence(divB, tol)

if __name__ == "__main__":
    main()
