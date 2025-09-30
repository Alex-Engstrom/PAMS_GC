# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 13:53:44 2025

@author: aengstrom
"""

import matplotlib.pyplot as plt
def make_pie_chart(slices: list, labels: list):
    fig1 = plt.figure(figsize = (5,5
                                 ))
    patches, texts, autotexts = plt.pie(slices, autopct='%1.1f%%', startangle=90)
    plt.axis('equal')
    plt.legend(patches, labels, loc="upper left")
    plt.show()

if __name__ == "__main__":
    labels = ['Valid Ambient Hours', 'QC Hours', 'Nulled Hours', 'Non-QC, Non-Ambient Hours']
    slices = [665, 77, 2,0]
    make_pie_chart(slices = slices, labels = labels)
        