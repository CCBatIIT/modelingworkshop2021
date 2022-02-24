---
layout: default
title: About
---

The workshop will be a one-week event comprising 10 periods (~3 hours each) that include both theory and computer laboratory. It is based on a [workshop in 2021](https://ccbatiit.github.io/modelingworkshop/index2021.html). The workshop will be primarily taught by David Minh. Some modules will be taught by guest lecturers.

Topics are as follows:

| Day | Period 1 | Period 2 |
| --- | -------- | -------- |
{% for day in (1..5) %} | {{ day }} | {% for period in (1..2) %} {% for module in site.data.modules %} {% if module.day == day %} {% if module.period == period %}{%- if module.bold %}<b>{% endif %}{{ module.title }}{% if module.bold %}</b>{% endif -%}{% if module.teacher %} ({{ module.teacher }}){% endif -%}{% if module.description %}. <i>{{ module.description }}</i> {% endif -%}
{% if module.basename != "None" %} [{% if module.slides == "ppt" %}[ppt](https://github.com/CCBatIIT/modelingworkshop/raw/main/slides/2022/{{ module.basename }}.ppt)/{% elsif module.slides == "pdf"%}{% else %}[key](https://github.com/CCBatIIT/modelingworkshop/raw/main/slides/2022/{{ module.basename }}.key)/{% endif -%}[pdf](https://github.com/CCBatIIT/modelingworkshop/raw/main/slides/2022/{{ module.basename }}.pdf)]. {% else %}.
{%- endif %} {%- endif %} {% endif %} {% endfor %} | {% endfor %}
{% endfor %}

Lectures may be recorded and computer laboratory exercises may be posted online for participants to complete asynchronously and as a free online resource.

We will soon update a list of software and web servers and prepare a virtual machine to download.


# Who?

This workshop will introduce advanced undergraduate students, graduate students, and research scientists in chemistry, biology, and related fields to computational methods for modeling biological macromolecules. Workshop participants will mostly be from Colombia, Romania, and the United States of America.

Up to 10 students may participate at Simon Bolivar University in Barranquilla, Colombia, and additional 10 students may participate online and receive support from the instructors.

Participants should have a laptop computer capable of running the virtual machine.

If you would like to be a formal participant in the workshop, please submit an [application](https://forms.gle/KgUUzhWsaF6XkkJ49). Formal participants will be able to attend the workshop live and participate in live discussions. Priority for formal participation will be given to affiliates of Simon Bolivar University and Illinois Tech and individuals working with a workshop instructor or collaborator.

After the workshop, lecture recordings and slides may be made freely available online.

# When?

The workshop will be held from March 14 to March 18, 2022. Modules will start at the following times, which are in Colombia Time (COT)/Central Daylight Time (CDT) and Eastern European Time (EST):

| Module | CDT  | EST  |
| ------ | ---- | ---- |
| 1      | 9 am | 4 pm |
| 2      | 2 pm | 9 pm |

All participants are welcome to work through the modules synchronously. Participants in COT/CDT are especially encouraged to participate live. However, it will probably make sense for remote participants, based on their time zone, to complete modules at a more convenient time. All participants can reach out to all workshop instructors for online office hours.

<!--
Online office hours will be held at the following times:

Workshop instructors will hold online office hours at the following times:

| Instructor | CDT  | EEST | ICT  |
| ---------- | ---- | ---- | ---- |
| David      | 7 pm | 3 am (+1 d) | 7 am (+1 d) |
| Laurentiu  | 1 am | 9 am | 1 pm |
| Soohaeng   | 8 am | 4 pm | 8 pm |
-->

# Where?

The workshop will be held at Simon Bolivar University in Barranqulla, Colombia.

# Why?

The purpose of this workshop will be
- to further the education of workshop participants in the field of computational biochemistry;
- to promote international scientific collaboration, especially between Colombia, Romania, and the United States; and
- to provide a rich cultural exchange experience that will foster this collaboration.
