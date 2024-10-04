Quantenteilchen verhalten sich scheinbar sehr unberechenbar. Bei einem Ballwurf ist
es einfach, die Flugbahn oder Auftreffort zu berechnen, nicht so in der Quantentheorie!
Denn ein Quantenteilchen kann alle erdenklichen Wege benutzen um von einem Ort zum
anderen zu kommen, und man kann die genaue Bewegung nicht so einfach nachvollziehen. Aber wie ist es möglich trotzdem eine Angabe darüber machen, welchen Weg es
gegangen ist um vom Abwurfpunkt bis zum Zielort zu kommen? Wo wird das Teilchen
am wahrscheinlichsten auftreffen? Für eine Antwort auf diese und weitere Fragen zum
Verhalten der Teilchen, gebraucht man das Pfadintegral, maßgeblich entwickelt durch
Richard Feynman. Dieses beschreibt genau die Summe über alle möglichen Pfade zwischen zwei Punkten und hilft dabei einige Rechnungen der Quantenmechanik klarer
und systematischer zu gestalten.Feynmans Formulierung des quantenmechanischen Pfadintegrals 
ist analytisch schwer bis garnicht lösbar. Wenn man sie jedoch mittels der Wick-Rotation 
vom Komplexen ins Reelle transformiert, lassen sich schnell Parallelen zur statistischen 
Mechanik erkennen. Das hat Physiker zu der Idee verholfen sich die Monte Carlo Methoden 
zu Nutze zu machen, was eine auf Wahrscheinlichkeiten basierende numerische Lösung des 
Pfadintegrals bedeutet.

Monte Carlo Methoden beschreiben eine breite Klasse von Algorithmen beschreiben, bei denen 
mit Hilfe von Wahrscheinlichkeitstheorie, das Pfadintegral numerisch gelöst werden kann. Grundlegende Abschnitte
sind dabei die Auffindung eines geeignete Wahrscheinlichkeitsmodells, mit dessen Verteilungsgesetzen eine Kette von Zufallszahlen erzeugt wird, die dann weiter genutzt werden
kann um Schätzwerte des Ausgangsproblems zu ermitteln. Zu dieser Klasse von Lösungsmethoden
gehört auch der sogenannte Metropolis Algorithmus, eingeführt durch Nicholas Metropolis. 
Er basiert auf importance sampling, wobei Punkte eines Phasenraums nicht komplett zufällig ausgewählt,
sondern vermehrt aus Regionen mit Zuständen, die entscheidende Beiträge zum Integral
liefern. Das heißt im Falle des Pfadintegrals, Pfade mittels der Wahrscheinlichkeitsverteilung die über das Pfadintegral hergeleitet wird, auszuwählen. Möglich wird das indem
man einen sogenannten Markov Prozess benutzt um die entscheidenden Zustände zu generieren. Dieser erzeugt eine Markov Kette, eine Kette möglicher Ereignisse, wobei jedes
zukünftige Kettenglied nur vom gegenwärtigen abhängt.
Erzeugt werden die Ereignisse durch einen Schritt im Metropolis Algorithmus.

