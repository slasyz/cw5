package main

import (
    "bufio"
    "fmt"
    "math"
    "os"
    "strconv"
    "strings"
)

const (
    eps   = 0.00000001
    spoil = 0.0000000001
    N     = 10
)

var (
    n int
)

func residual(x [N]float64, a [N][N]float64) (R [N]float64) {
    for i := 1; i <= n; i++ {
        R[i] = a[i][n+1]

        for j := 1; j <= n; j++ {
            R[i] -= a[i][j] * x[j]
        }
    }

    return
}

func jordan_one(a [N][N]float64) (x [N]float64) {
    /*for k := 1; k <= n; k++ {
          tmp := a[k][k]
          if math.Abs(float64(tmp)) < eps {
              fmt.Printf("Слишком маленькое значение a[%d][%d][%d]\n", k-1, k, k)
          }
          for j := k; j <= n+1; j++ {
              a[k][j] /= tmp
          }

          for i := k + 1; i <= n; i++ {
              tmp = a[i][k]
              for j := k; j <= n+1; j++ {
                  a[i][j] -= a[k][j] * tmp
              }
          }
      }

      for i := n; i >= 1; i-- {
          x[i] = a[i][n+1]
          for j := i + 1; j <= n; j++ {
              x[i] -= a[i][j] * x[j]
          }
      }*/

    for k := 1; k <= n; k++ {
        tmp := a[k][k]
        if math.Abs(tmp) < eps {
            fmt.Printf("Слишком маленькое значение a[%d][%d][%d] = %.16f\n", k-1, k, k, tmp)
        }
        for j := k; j <= n+1; j++ {
            a[k][j] /= tmp
        }

        for i := 1; i <= n; i++ {
            if i == k {
                continue
            }

            tmp = a[i][k]
            for j := k; j <= n+1; j++ {
                a[i][j] -= a[k][j] * tmp
            }
        }
    }

    /*for i := 1; i <= n; i++ {
        for j := 1; j <= n+1; j++ {
            fmt.Printf("%9.6f ", a[i][j])
        }
        fmt.Println()
    }*/

    for i := 1; i <= n; i++ {
        x[i] = a[i][n+1]
    }
    return
}

func gauss_select(a [N][N]float64) (x [N]float64) {
    for k := 1; k <= n; k++ {
        max := a[k][k]
        i_max := k

        // Ищем максимальный элемент в столбце, запоминаем его
        for i := k + 1; i <= n; i++ {
            if math.Abs(a[i][k]) > math.Abs(max) {
                max = a[i][k]
                i_max = i
            }
        }

        // Меняем i_max строку со строкой k
        if i_max != k {
            fmt.Printf("Меняем строки %d <-> %d\n", i_max, k)
            tmp := a[i_max]
            a[i_max] = a[k]
            a[k] = tmp
        }

        // Решаем всё как раньше
        tmp := a[k][k]
        if math.Abs(tmp) < eps {
            fmt.Printf("Слишком маленькое значение a[%d][%d][%d] = %.16f\n", k-1, k, k, tmp)
        }
        for j := k + 1; j <= n+1; j++ {
            a[k][j] /= tmp
        }

        for i := 1; i <= n; i++ {
            if i == k {
                continue
            }

            tmp = a[i][k]
            for j := k; j <= n+1; j++ {
                a[i][j] -= a[k][j] * tmp
            }
        }
    }

    for i := 1; i <= n; i++ {
        x[i] = a[i][n+1]
    }
    return
}

func inv(a [N][N]float64) (A_inv [N][N]float64) {
    for i := 1; i <= n; i++ {
        for j := 1; j <= n; j++ {
            if i == j {
                a[i][n+j] = 1
            } else {
                a[i][n+j] = 0
            }
        }
    }

    for k := 1; k <= n; k++ {
        tmp := a[k][k]
        if math.Abs(tmp) < eps {
            fmt.Printf("Слишком маленькое значение a[%d][%d][%d] = %.16f\n", k-1, k, k, tmp)
        }
        for j := k + 1; j <= 2*n; j++ {
            a[k][j] /= tmp
        }

        for i := 1; i <= n; i++ {
            if i == k {
                continue
            }

            tmp = a[i][k]
            for j := k; j <= 2*n; j++ {
                a[i][j] -= a[k][j] * tmp
            }
        }
    }

    return a
}

func main() {
    scanner := bufio.NewScanner(os.Stdin)
    fmt.Print("Введите n: ")
    scanner.Scan()
    n, _ = strconv.Atoi(scanner.Text())

    var A [N][N]float64
    fmt.Println("Введите матрицу A: ")
    for i := 1; i <= n; i++ {
        scanner.Scan()
        line_str := strings.Split(strings.Trim(scanner.Text(), " "), " ")
        for j := 1; j <= n; j++ {
            tmp_float, _ := strconv.ParseFloat(line_str[j-1], 32)
            A[i][j] = float64(tmp_float)
        }
    }

    fmt.Println("Введите столбец b: ")
    scanner.Scan()
    line_str := strings.Split(strings.Trim(scanner.Text(), " "), " ")
    for i := 1; i <= n; i++ {
        tmp_float, _ := strconv.ParseFloat(line_str[i-1], 32)
        A[i][n+1] = float64(tmp_float)
    }

    C := A
    C[1][1] = A[1][1] * spoil

    fmt.Println()
    fmt.Println("Воспользуемся методом единственного деления: ")
    x := jordan_one(A)
    fmt.Printf("x = %19.16f\n", x[1:n+1])
    R := residual(x, A)
    fmt.Printf("R = %19.16f\n", R[1:n+1])

    fmt.Println()
    fmt.Println("Испортим немного матрицу и снова решим систему: ")
    x = jordan_one(C)
    fmt.Printf("x = %19.16f\n", x[1:n+1])
    R = residual(x, C)
    fmt.Printf("R = %19.16f\n", R[1:n+1])

    fmt.Println()
    fmt.Println("Воспользуемся методом с выбором главного элемента (по столбцу): ")
    x = gauss_select(A)
    fmt.Printf("x = %19.16f\n", x[1:n+1])
    R = residual(x, A)
    fmt.Printf("R = %19.16f\n", R[1:n+1])

    fmt.Println()
    fmt.Println("Испортим немного матрицу и снова решим систему: ")
    x = gauss_select(C)
    fmt.Printf("x = %19.16f\n", x[1:n+1])
    R = residual(x, C)
    fmt.Printf("R = %19.16f\n", R[1:n+1])

    fmt.Println()
    fmt.Println("Найдём обратную к матрице A: ")
    A_inv := inv(A)
    for i := 1; i <= n; i++ {
        for j := n + 1; j <= 2*n; j++ {
            fmt.Printf("%9.6f ", A_inv[i][j])
        }
        fmt.Println()
    }
}
