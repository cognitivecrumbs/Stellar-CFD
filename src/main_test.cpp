// // // #include <iostream>
// // // #include <vector>

// // // class MyClass {
// // // public:
// // //     int value;

// // //     MyClass(int val) : value(val) {
// // //         std::cout << "MyClass constructor: " << value << std::endl;
// // //     }

// // //     ~MyClass() {
// // //         std::cout << "MyClass destructor: " << value << std::endl;
// // //     }

// // //     void display() const {
// // //         std::cout << "Value: " << value << std::endl;
// // //     }
// // // };

// // // int main() {
// // //     // Step 1: Create a vector to store raw pointers to MyClass objects
// // //     std::vector<MyClass*> objectVector;

// // //     // Step 2: Add objects to the vector (using new)
// // //     objectVector.push_back(new MyClass(10));
// // //     objectVector.push_back(new MyClass(20));
// // //     objectVector.push_back(new MyClass(30));

// // //     // Step 3: Access and display the objects
// // //     std::cout << "\nObjects in the vector:" << std::endl;
// // //     for (const auto& obj : objectVector) { 
// // //         obj->display();
// // //     }

// // //     // Step 4: Clean up memory manually
// // //     std::cout << "\nCleaning up memory:" << std::endl;
// // //     for (auto& obj : objectVector) {
// // //         delete obj;  // Free memory for each object
// // //         obj = nullptr; // Avoid dangling pointers
// // //     }

// // //     // Optional: Clear the vector (since it only contains pointers)
// // //     objectVector.clear();

// // //     return 0;
// // // }


// // #include <iostream>
// // #include <vector>

// // class MyClass {
// // public:
// //     int value;

// //     MyClass(int val) : value(val) {
// //         std::cout << "MyClass constructor: " << value << std::endl;
// //     }

// //     ~MyClass() {
// //         std::cout << "MyClass destructor: " << value << std::endl;
// //     }

// //     void display() const {
// //         std::cout << "Value: " << value << std::endl;
// //     }
// // };

// // int main() {
// //     const int SIZE = 3;

// //     // Step 1: Preallocate a vector with nullptrs
// //     std::vector<MyClass*> objectVector(SIZE, nullptr);

// //     // Step 2: Assign actual objects into the preallocated slots
// //     for (size_t i = 0; i < objectVector.size(); ++i) {
// //         objectVector[i] = new MyClass((i + 1) * 10); // Assign new objects
// //     }

// //     // Step 3: Print the vector pointer addresses and the object values
// //     std::cout << "\nObjects in the vector (pointer addresses and values):" << std::endl;
// //     for (size_t i = 0; i < objectVector.size(); ++i) {
// //         std::cout << "Pointer Address: " << objectVector[i] << " -> ";
// //         objectVector[i]->display(); // Print object value
// //     }

// //     // Step 4: Clean up memory manually
// //     std::cout << "\nCleaning up memory:" << std::endl;
// //     for (auto& obj : objectVector) {
// //         delete obj;  // Free memory for each object
// //         obj = nullptr; // Avoid dangling pointers
// //     }

// //     // Optional: Clear the vector
// //     objectVector.clear();

// //     return 0;
// // }

// #include <iostream>
// #include <cmath>
// #include <chrono>
// #include <glm/glm.hpp>

// int main() {
//     const int N = 1000000; // Number of iterations
//     double result = 0.0;

//     // Measure sqrt(i)^3
//     auto start1 = std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec1(i,i);
//         float dist1 = std::pow(vec1.x,2)+std::pow(vec1.y,2);
//         // result += std::pow(std::sqrt(dist1), 3);
//         dist1 = std::sqrt(dist1);
//         result += dist1*dist1*dist1;
//     }
//     auto end1 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for sqrt(i)^3: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count() 
//               << " microseconds\n";

//     result = 0.0;

//     // Measure i*sqrt(i)
//     auto start2 = std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec2(i,i);
//         float dist2 = std::pow(vec2.x,2)+std::pow(vec2.y,2);
//         result += dist2 * std::sqrt(dist2);
//     }
//     auto end2 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for i*sqrt(i): " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count() 
//               << " microseconds\n";


//     // Measure i*sqrt(i)
//     auto start3 = std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec3(i,i);
//         float dist3 = std::sqrt(std::pow(vec3.x,2)+std::pow(vec3.y,2));
//         // result += dist2 * std::sqrt(dist2);
//     }
//     auto end3 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for squared: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3).count() 
//               << " microseconds\n";


//     // Measure i*sqrt(i)
//     auto start4= std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec4(i,i);
//         // float dist4 = vec4.length();
//         float dist4 = dot(vec4,vec4);
//         // result += dist2 * std::sqrt(dist2);
//     }
//     auto end4 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for len: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4).count() 
//               << " microseconds\n";

//     // Measure i*sqrt(i)
//     auto start5 = std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec5(i,i);
//         float dist5 = std::pow(vec5.x,2)+std::pow(vec5.y,2);
//         float dist53 = dist5*std::sqrt(dist5);
//         // result += dist2 * std::sqrt(dist2);
//     }
//     auto end5 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for squared: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end5 - start5).count() 
//               << " microseconds\n";


//     // Measure i*sqrt(i)
//     auto start6= std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec6(i,i);
//         // float dist6 = vec6.length();
//         float dist6 = sqrt(dot(vec6,vec6));
//         float dist63 = std::pow(dist6,3);
//         // result += dist2 * std::sqrt(dist2);
//     }
//     auto end6 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for len: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end6 - start6).count() 
//               << " microseconds\n";

//     // Measure i*sqrt(i)
//     auto start7= std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec7(i,i);
//         // float dist7 = vec7.length();
//         float dist7 = std::sqrt(dot(vec7,vec7));
//         float dist73 = dist7*dist7*dist7;
//         // result += dist2 * std::sqrt(dist2);
//     }
//     auto end7 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for len: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end7 - start7).count() 
//               << " microseconds\n";

//     // Measure i*sqrt(i)
//     auto start8= std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec8(i,i);
//         // float dist7 = vec7.length();
//         float dist8_s = dot(vec8,vec8);
//         float dist83 = dist8_s*std::sqrt(dist8_s);
//         // result += dist2 * std::sqrt(dist2);
//     }
//     auto end8 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for len: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end8 - start8).count() 
//               << " microseconds\n";

//     // Measure i*sqrt(i)
//     auto start9= std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec9(i,i);
//         // float dist7 = vec7.length();
//         float dist9_s = dot(vec9,vec9);
//         float dist93 = dist9_s*std::sqrt(dist9_s);
//         // result += dist2 * std::sqrt(dist2);
//     }
//     auto end9 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for len: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end9 - start9).count() 
//               << " microseconds\n";

//     // Measure i*sqrt(i)
//     auto start8= std::chrono::high_resolution_clock::now();
//     for (int i = 1; i <= N; ++i) {
//         glm::vec2 vec8(i,i);
//         // float dist7 = vec7.length();
//         float dist8_s = dot(vec8,vec8);
//         float dist83 = dist8_s*std::sqrt(dist8_s);
//         // result += dist2 * std::sqrt(dist2);
//     }
//     auto end8 = std::chrono::high_resolution_clock::now();
//     std::cout << "Time for len: " 
//               << std::chrono::duration_cast<std::chrono::microseconds>(end8 - start8).count() 
//               << " microseconds\n";
//     return 0;
// }


#include <iostream>
#include <cmath>
#include <chrono>

int main() {
    double value = 2.5;
    int iterations = 100000000;

    // Time sqrt()
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        double result = std::sqrt(value);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto sqrtDuration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Time value * value
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        double result = value * value;
    }
    end = std::chrono::high_resolution_clock::now();
    auto multDuration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "sqrt(): " << sqrtDuration.count() << " microseconds\n";
    std::cout << "value * value: " << multDuration.count() << " microseconds\n";

    return 0;
}