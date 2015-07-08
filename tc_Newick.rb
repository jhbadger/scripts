#!/usr/bin/env ruby

require 'test/unit'
require 'Newick'

class TC_Newick < Test::Unit::TestCase
  
  def setup
    @tree = NewickTree.new("((((Cat:0.5,Dog:0.5):2.0,Bear:2.0):3.0,Ape:3.0):10.0,Earthworm:4.0);")
  end

  def test_to_s
    assert_equal("((((Cat:0.5,Dog:0.5):2.0,Bear:2.0):3.0,Ape:3.0):10.0,Earthworm:4.0);", @tree.to_s, "to_s error")
  end

  def test_reroot
    @tree.reroot(@tree.findNode("Dog"))
    assert_equal("(Dog:0.25,(Cat:0.5,(Bear:2.0,(Earthworm:14.0,Ape:3.0):3.0):2.0):0.25);", @tree.to_s, "reroot error")
    @tree.reroot(@tree.findNode("Ape"))
    assert_equal("(Ape:1.5,(Earthworm:14.0,(Bear:2.0,(Dog:0.5,Cat:0.5):2.0):3.0):1.5);", @tree.to_s, "repeat reroot error")
  end

  def test_unroot
    @tree.unroot
    assert_equal("(Earthworm:14.0,((Cat:0.5,Dog:0.5):2.0,Bear:2.0):3.0,Ape:3.0);", @tree.to_s, "unroot error")
  end

  def test_midroot
    @tree.midpointRoot
    assert_equal("(Earthworm:4.25,(((Cat:0.5,Dog:0.5):2.0,Bear:2.0):3.0,Ape:3.0):9.75);", @tree.to_s, "midpoint error")
  end

  def test_lca
    wormNode = @tree.findNode("Earthworm")
    bearNode = @tree.findNode("Bear")
    lca = wormNode.lca(bearNode)
    assert_equal(4.0, wormNode.distToAncestor(lca), "worm lca dist error")
    assert_equal(15.0, bearNode.distToAncestor(lca), "bear lca dist error")
  end
 
  def test_dMatrix
    distMatrix = @tree.distanceMatrix
    assert_equal(8.5, distMatrix["Cat"]["Ape"], "dist matrix error")
  end

  def test_clades
    clades = @tree.clades
    assert_equal(3, clades.size, "clades error")
  end

  def teardown
  end

end
