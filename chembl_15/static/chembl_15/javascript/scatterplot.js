

        // Width and height
        var w =300;
        var h = w;
        var pad = {left:90, top:40, right:50, bottom:70}
        var innerpad = 5
        
        var totalh = h + pad.top + pad.bottom
        var totalw = w + pad.left + pad.right

         
        var svg = d3.select("body")
                    .append("svg")
                    .attr("width", totalw)
                    .attr("height", totalh);

        // Panel for dot plot.
        var dotplot = svg.append("g")
                    .attr("id", "dotplot")
                    .attr("transform", "translate(" + pad.left +  "," + pad.top +  ")")
                    
        // Scales

        var xScale = d3.scale.linear()
            //.domain([d3.min(top_acts, function(d){return d[1];}), d3.max(top_acts, function(d){return d[1];})])
            .domain([3,10])
            .range([20,w-20])
            .clamp(true);

        var yScale = d3.scale.ordinal()
            .domain(top_mols)
            .rangeBands([20,h], .1);

        // Background
        dotplot.append("rect")
             .attr("height", h)
             .attr("width", w)
             .attr("fill", d3.rgb(200, 200, 200))
             .attr("stroke", "black")
             .attr("stroke-width", 1)
             .attr("pointer-events", "none");


        // Axis labels
        dotplot.append("text")
               .attr("id", "xaxis")
               .attr("class", "axes")
               .attr("x", w/2)
               .attr("y", h+pad.bottom*0.7)
               .text("pIC50")
               .attr("dominant-baseline", "middle")
               .attr("text-anchor", "middle")
               .attr("fill", "slateblue");
        dotplot.append("text")
               .attr("id", "yaxis")
               .attr("class", "axes")
               .attr("x", -pad.left*0.8)
               .attr("y", h/2)
               .text("molregno")
               .attr("dominant-baseline", "middle")
               .attr("text-anchor", "middle")
               .attr("transform","rotate(270, " + -pad.left*0.9 + "," + h/2 +")")
               .attr("fill", "slateblue");

        // Axis scales
        var xticks = xScale.ticks(5);
        var yticks = top_mols;

        dotplot.selectAll("empty")
               .data(xticks)
               .enter()
               .append("text")
               .attr("class", "axes")
               .text(function(d){
                            return d3.format(".2f")(d);})
               .attr("x", function(d){
                                 return xScale(d);
                                 })
               .attr("y", h+pad.bottom*0.3)
               .attr("dominant-baseline", "middle")
               .attr("text-anchor", "middle");

        dotplot.selectAll("empty")
               .data(yticks)
               .enter()
               .append("text")
               .attr("class", "axes")
               .text(function(d){
                                return d;
                                })
               .attr("x", -pad.left*0.05)
               .attr("y", function(d){ 
                                    return yScale(d);
                                    })
               .attr("dominant-baseline", "middle")
               .attr("text-anchor", "end");

        dotplot.selectAll("empty")
               .data(xticks)
               .enter()
               .append("line")
               .attr("class", "axes")
               .attr("x1", function(d){
                                    return xScale(d);
                                    })
               .attr("x2", function(d){
                                    return xScale(d);
                                    })
               .attr("y1", 0)
               .attr("y2", h)
               .attr("stroke", "white")
               .attr("stroke-width", 1);

        dotplot.selectAll("empty")
               .data(yticks)
               .enter()
               .append("line")
               .attr("class", "axes")
               .attr("y1", function(d){
                                    return yScale(d);
                                    })
               .attr("y2", function(d){
                                    return yScale(d);
                                    })
               .attr("x1", 0)
               .attr("x2", w)
               .attr("stroke", "white")
               .attr("stroke-width", 1);
        // Add dots.  
        dotplot.selectAll("empty")
               .data(top_acts)
               .enter()
               .append("circle")
               .attr("cy", function(d) {
                            return yScale(d[0])+ Math.floor(Math.random()*13) - 6;
                            })
               .attr("cx", function(d) {
                                return xScale(d[1]);
                                })
               .attr("r", 3)
               .attr("stroke", "black")
               .attr("stroke-width", 1)
               .attr("fill", "red")
               .on("mouseover", function(d){
                   d3.select(this).attr("stroke", "white")
                   dotplot.append("foreignObject").attr("id", "prefname")
                         .attr("width","150")
                         .attr("height","300")
                         .attr("x",xScale(d[1]))
                         .attr("y", yScale(d[0]))
                         .append("xhtml:body")
                         .append("p")
                         .text(d[4])
                         .attr("dominant-baseline", "middle")
                         .attr("text-anchor", "middle")
                    dotplot.append("image").attr("id","2dgraph")
                         .attr("xlink:href", "/static/chembl_15/png/test.png")
                         .attr("width", 50)
                         .attr("height", 50)
                         .attr("x",xScale(d[1]) + 1 * 30)
                         .attr("y", yScale(d[0]) - (1 + 0.35) * 20) 
                    })
                .on("mouseout", function(d){
                    d3.selectAll("foreignObject#prefname").remove()
                    d3.selectAll("image").remove()
                    d3.select(this).attr("stroke","black")
                    });
