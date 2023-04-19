# tf
Transfer Function for Linear System

### 순서
1. cpx.py
2. tf.py
   - Class Tfalg()
   - Class Tf()


***

## 1. cpx.py
> This is for calculation of Complex Number (Sum, Product (with list))
<br> 복소수 계산을 위한 클래스입니다. (합, 곱 (with 여러개는 리스트로 받음))</br>

- comp_sum(a,b): Sum of Complex Number
- comp_sum_ls(list): (Multiple) Sum of Complex Number
- comp_mul(a,b): Product of Complex Number</li>
- comp_mul_ls(list): (Multiple) Product of Complex Number</li>


## 2. tf.py
### Class Tfalg()
> This is just a Polynomials with Rational Coefficients
<br> 유리수계수 다항함수입니다. </br>

#### Class Variables
  - root_dict
  <br> 다항함수의 근을 가진 Dictionary (근 : 중복횟수) </br>
  - n_coef
  <br> 최고차항 계수 </br>
  - order
  <br> 다항함수의 차수 </br>
  - coef
  <br> 다항함수의 계수(내림차순) </br>
  
#### Class Functions
  - setCoef(coef)
  <br> 다항함수의 계수를 set **(이 때, 근을 자동으로 구해주지 않습니다!!!)**
  <br> 그러므로, 근은 아래의 setRoots를 이용해서 따로 세팅합니다. </br>
  - setRoots(roots,n)
  <br> 다항함수의 근을 리스트로 받아서 근을 세팅합니다. (n은 최고차항 계수) **(이 때, 계수는 calc_coef를 통해 자동으로 구해집니다!!)*** </br>
  - getRoots()
  <br> 다항함수의 근을 반환합니다. </br>
  - getCoef()
  <br> 다항함수의 계수를 반환합니다. </br>
  - calc_coef()
  <br> 갖고 있는 근을 가지고 계수를 계산합니다. </br>
  - calc(x)
  <br> 주어진 다항함수에 x를 넣었을 때, 얼마가 나오는지 계산합니다. </br>
  - find_roots(step)
  <br> 계수가 주어졌을 때, 다항함수의 근을 계산해줍니다. (단, 수치적 계산으로 오차가 있을 수 있습니다. Step은 반복 횟수입니다.) </br>
  <br></br>
     Gradient Descent 방법을 이용해서 다항함수의 근을 계산합니다.
        
       ![linsys_code1](https://user-images.githubusercontent.com/40926406/233068014-c9a7fcb1-5a26-4f29-89ea-65eb15441cc9.png)
       ![linsys_code2](https://user-images.githubusercontent.com/40926406/233068036-cf8b7cae-a68d-4b9a-8a18-f7c34dfcb7b9.png)
        
     이 방법을 이용하여 근 하나를 찾고, 하나를 찾으면 차수를 내려서 다시 근을 찾는 것을 반복합니다.
     <br> 1차식, 2차식일 때는 근의 공식을 이용해서 근을 찾습니다.</br>
     
  - tfalg_print()
  <br> 주어진 다항함수를 print합니다. </br>
  
### Class Tf()
> This is a Transfer Function with Linear System
<br> 선형시스템의 전달함수입니다. (복소수함수의 분수함수꼴) </br>

#### Class Variables
  - num
  <br> 전달함수의 분자부분 </br>
  - den
  <br> 전달함수의 분모부분 </br>
  
#### Class Functions
  - setZeros(zeros), setPoles(poles)
  <br> 전달함수의 zero와 pole을 set합니다. 이 때, pole,zero에 겹치는 해를 자동으로 제거합니다. </br>
  - getZeros(), getPoles()
  <br> 전달함수의 zero와 pole을 반환합니다.</br>
  - elim()
  <br> 전달함수 분자 분모에 겹치는 해를 제거하기 위해서 만든 함수입니다.. </br>
  - tf_divide(), tf_2nd_divide()
  <br> 전달함수를 부분분수 꼴로 다 쪼개서 tf의 리스트를 반환합니다.
  <br> tf_2nd_divide() 는 2차식의 분모를 가진 전달함수를 라플라스 변환하기 쉽도록 cos파트와 sin 파트로 나눕니다.(계수만 계산합니다) </br>
  
    ![linsys_code3](https://user-images.githubusercontent.com/40926406/233073209-931efedc-8a08-4884-80ab-e0b7c859f0e6.png)
  
     <br></br>
     이 때, 각 전달함수에 들어가는 분자 부분은 다음과 같이 계산합니다.
  
    ![linsys_code4](https://user-images.githubusercontent.com/40926406/233075443-83a09173-08e2-413e-b0df-0118fbcdb491.png)
    
    ![linsys_code5](https://user-images.githubusercontent.com/40926406/233075462-7b3a862d-b21f-43b0-b7b9-b47841f2a0cb.png)

  - tf_response()
  <br> 라플라스 변환을 이용해서 time-domain의 Response를 반환합니다. **(이 때, tf에 적절한 Input도 넣어서 계산해야 합니다!)** </br>

